# Import necessary libraries
import json
import scanpy as sc
import os
import sys
import subprocess
import time
import numpy as np
import pandas as pd
from multiprocessing import Process
from pathlib import Path

import cnmf
from cnmf import cNMF


def validate_config(config):
    """
    Checks the configuration for potential issues. Exits if critical errors are found.
    This function validates the parameters from the JSON config file.
    """
    print("--- Validating Configuration ---")
    errors = []

    # Check 1: Ensure all required keys are present in the JSON file.
    required_keys = [
        "output_dir", "output_name", "ad_file", "component_list",
        "num_iter", "total_nodes", "threads_per_node"
    ]
    for key in required_keys:
        if key not in config:
            errors.append(f"Missing required key in config file: '{key}'")
    
    if any("Missing required key" in e for e in errors):
        for error in errors:
            print(f"ERROR: {error}")
        print("--------------------------------\n")
        sys.exit(1)

    # Check 2: Validate file and directory paths.
    ad_file_path = Path(config["ad_file"])
    if not ad_file_path.exists():
        errors.append(f"Input data file/directory not found: {ad_file_path}")
    else:
        # Check if the path is a directory (10x MTX), .h5ad, or .h5 (10x HDF5).
        if not (ad_file_path.is_dir() or ad_file_path.suffix in ['.h5ad', '.h5']):
            errors.append(f"Input data '{ad_file_path}' is not a recognized format (directory, .h5ad, or .h5).")

    # Check 3: Check output directory and create it if it doesn't exist.
    output_dir_path = Path(config["output_dir"])
    try:
        output_dir_path.mkdir(parents=True, exist_ok=True)
        print(f"Output directory checked/created: {output_dir_path}")
    except OSError as e:
        errors.append(f"Could not create output directory {output_dir_path}: {e}")

    # Check 4: Validate numeric parameters.
    if not (isinstance(config["num_iter"], int) and config["num_iter"] > 0):
        errors.append(f"'num_iter' must be a positive integer. Got: {config['num_iter']}")
    if not (isinstance(config["total_nodes"], int) and config["total_nodes"] > 0):
        errors.append(f"'total_nodes' must be a positive integer. Got: {config['total_nodes']}")
    if not (isinstance(config["threads_per_node"], int) and config["threads_per_node"] > 0):
        errors.append(f"'threads_per_node' must be a positive integer. Got: {config['threads_per_node']}")

    # Check 5: Validate the list of components.
    if not isinstance(config["component_list"], list) or not all(isinstance(k, int) and k > 0 for k in config["component_list"]):
        errors.append(f"'component_list' must be a list of positive integers. Got: {config['component_list']}")

    if errors:
        print("\nConfiguration validation failed with the following errors:")
        for error in errors:
            print(f"- {error}")
        print("--------------------------------\n")
        sys.exit(1)
    else:
        print("Configuration seems OK.")
        print("--------------------------------\n")


def load_anndata_for_check(path_str):
    """
    Loads AnnData object from various formats into memory.
    """
    data_path = Path(path_str)
    print(f"--- Loading input data from {data_path} ---")

    adata = None
    if data_path.is_dir():
        print("Detected directory, loading as 10x MTX format.")
        try:
            adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', cache=False)
            print("Successfully loaded 10x MTX data.")
        except Exception as e:
            print(f"ERROR: Failed to load 10x MTX data: {e}")
            sys.exit(1)
    elif data_path.is_file() and data_path.suffix == '.h5ad':
        print("Detected .h5ad file, loading.")
        try:
            adata = sc.read_h5ad(data_path)
            print("Successfully loaded h5ad data.")
        except Exception as e:
            print(f"ERROR: Failed to load h5ad file: {e}")
            sys.exit(1)
    elif data_path.is_file() and data_path.suffix == '.h5':
        print("Detected .h5 file, loading as 10x HDF5 format.")
        try:
            adata = sc.read_10x_h5(data_path)
            print("Successfully loaded 10x HDF5 data.")
        except Exception as e:
            print(f"ERROR: Failed to load 10x HDF5 file: {e}")
            sys.exit(1)
    else:
        print(f"ERROR: Unrecognized input data format for {data_path}.")
        sys.exit(1)
    
    return adata

def run_finalize(cnmf_object, num_k):
    """
    Function to combine and build consensus matrices for a given K.
    """
    print(f"Starting finalization for K={num_k}...")
    cnmf_object.combine(components=num_k)
    cnmf_object.consensus(k=num_k, density_threshold=0.2)
    print(f"Finalization for K={num_k} completed.")

def main():
    """
    Main execution logic of the script.
    """
    if len(sys.argv) < 2:
        print("Usage: python your_script_name.py <path_to_config.json>")
        sys.exit(1)
    
    config_file = sys.argv[1]
    if not Path(config_file).is_file():
        print(f"ERROR: Config file not found at {config_file}")
        sys.exit(1)

    with open(config_file, "r") as json_file:
        config = json.load(json_file)
    
    print("--- Configuration Parameters ---")
    for key, value in config.items():
        print(f"{key}: {value}")
    print("------------------------------\n")
    
    validate_config(config)

    output_dir = config["output_dir"]
    output_name = config["output_name"]
    ad_file = config["ad_file"]
    
    component_list = config["component_list"]
    num_genes = config.get("num_genes",2000)
    num_iter = config["num_iter"]
    total_nodes = config["total_nodes"]
    threads_per_node = config["threads_per_node"]
    
    script_dir = Path(__file__).parent
    factor_script_file = config["factor_script_file"]
    print(f"factorize script: {factor_script_file}")
    
    # Detect cnmf.py path via Python import (most reliable)
    try:
        # Get the file path of the imported module
        cnmf_module_path = Path(cnmf.__file__).resolve()
    
        # If __file__ points to __init__.py, move to cnmf.py
        if cnmf_module_path.name == "__init__.py":
            candidate_path = cnmf_module_path.parent / "cnmf.py"
        else:
            candidate_path = cnmf_module_path
    
        if not candidate_path.exists():
            raise FileNotFoundError(f"cnmf.py not found: {candidate_path}")
    
        cNMF_python_path = str(candidate_path)
        print(f"cNMF_path: {cNMF_python_path}")
    
    except ImportError:
        print("ERROR: cnmf is not installed in the current Python environment.")
        sys.exit(1)
    
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    
    # Load the data into an AnnData object, regardless of original format.
    ad_data = load_anndata_for_check(ad_file)
    
    # Make variable names unique to prevent issues downstream.
    ad_data.var_names_make_unique()
    print("Made variable names unique.")

    print(f"Number of cells: {ad_data.shape[0]}")
    print(f"Number of genes: {ad_data.shape[1]}\n")
    
    if num_genes=="All" or num_genes=="all":
        #If "all" is specified, use all genes for calculation
        num_genes = ad_data.shape[1]
        print(f"Use all genes for program calculation. Num of genes:{num_genes}")
    if num_genes>ad_data.shape[1]:
        num_genes = ad_data.shape[1]
        print(f"Num of high var genes is more than dataset. Num of genes:{num_genes}")
        
    # Initialize the cNMF object. This also creates the main output directories.
    cnmf_obj = cNMF(output_dir=output_dir, name=output_name)
    
    # This will be the path we pass to cNMF. It defaults to the original file path.
    cnmf_input_path = ad_file
    
    # Check if the original input was a .h5 file.
    if Path(ad_file).suffix == '.h5':
        # Define a path for the temporary .h5ad file inside the output directory.
        temp_h5ad_path = os.path.join(output_dir,
                                      output_name,
                                      'cnmf_tmp',
                                      f"{output_name}.input_data.h5ad")
        
        print(f"Input is a .h5 file. Converting to a temporary .h5ad file for cNMF: {temp_h5ad_path}")
        # Write the AnnData object to the new .h5ad file.
        ad_data.write_h5ad(temp_h5ad_path)
        
        # Update the path to be used by cNMF.
        cnmf_input_path = temp_h5ad_path

    # We no longer need the in-memory AnnData object.
    del ad_data

    total_threads = total_nodes * threads_per_node
    
    print("--- Preparing cNMF analysis ---")
    # Pass the appropriate file path (original or temporary .h5ad) to the prepare function.
    cnmf_obj.prepare(counts_fn=cnmf_input_path,
                     components=component_list,
                     num_highvar_genes=num_genes,
                     n_iter=num_iter, seed=1)
    print("cNMF preparation complete.\n")

    # The rest of the script remains unchanged.
    print("--- Submitting factorization jobs ---")
    output_file_name = os.path.join(output_dir, "cnmf_factor_%j")
    sbatch_output_info = f"-o {output_file_name}.out -e {output_file_name}.err"

    for i in range(total_nodes):
        command_factor = (
            f"sbatch {sbatch_output_info} {factor_script_file} "
            f"{total_threads} {threads_per_node} {i * threads_per_node} "
            f"'{output_dir}' '{output_name}' '{cNMF_python_path}'"
        )
        print(f"Submitting job {i+1}/{total_nodes}: {command_factor}")
        subprocess.run(command_factor, shell=True)
    print("All jobs submitted.\n")

    res_folder = os.path.join(output_dir, output_name, "cnmf_tmp")
    
    print("--- Monitoring job progress ---")
    finish_flag = dict(zip(component_list, [False] * len(component_list)))
    start = time.time()
    
    while not all(finish_flag.values()):
        elapsed_time = np.round(time.time() - start, 2)
        print(f"[{elapsed_time}s] Checking progress...")
        
        for target_k in component_list:
            if finish_flag[target_k]:
                continue
            
            search_command = f"ls {res_folder}/{output_name}.*.k_{target_k}.* 2>/dev/null | wc -l"
            try:
                res_wc = subprocess.check_output(search_command, shell=True, text=True)
                num_line = int(res_wc.strip())
            except (subprocess.CalledProcessError, ValueError):
                num_line = 0

            if num_line >= num_iter:
                if not finish_flag[target_k]:
                    print(f"*** K={target_k}: All {num_iter} iterations finished. Starting finalization. ***")
                    background_process = Process(target=run_finalize, args=(cnmf_obj, target_k), daemon=True)
                    background_process.start()
                    finish_flag[target_k] = True
            else:
                print(f"K={target_k}: {num_line}/{num_iter} iterations completed.")
        
        if all(finish_flag.values()):
            break
        
        print()
        time.sleep(30)

    print("\n--- All factorization jobs completed. Waiting for finalization... ---")
    start_finalize_wait = time.time()
    while True:
        elapsed_time = np.round(time.time() - start_finalize_wait, 2)
        stat_folder = os.path.join(output_dir, output_name)
        
        search_command = f"ls {stat_folder}/{output_name}.gene_spectra_score.k_*.txt"+\
            f" 2>/dev/null | wc -l"
        try:
            res_wc = subprocess.check_output(search_command, shell=True, text=True)
            num_line = int(res_wc.strip())
        except (subprocess.CalledProcessError, ValueError):
            num_line = 0

        if num_line >= len(component_list):
            print("All K values have been processed. Finalization complete.")
            break
        else:
            print(f"[{elapsed_time}s] Combine & consensus: {num_line}/{len(component_list)} finished.")
        time.sleep(5)

    print("\n--- Generating K selection plot ---")
    try:
        cnmf_obj.k_selection_plot()
        print("Analysis complete. K selection plot has been generated.")
    except Exception as e:
        print(f"Could not generate K selection plot. Error: {e}")
        print("The analysis is otherwise complete. You can try generating the plot manually.")

if __name__ == "__main__":
    main()