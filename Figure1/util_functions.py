import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
import multiprocessing as mp
from itertools import combinations

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

import torch
import gc

import scanpy as sc
import pickle
import os
import torch

def load_files(input_file,sgRNA_file,dict_file,pca_file,obsm_key):
    print("read input")
    if os.path.isfile(pca_file):
        X = pd.read_pickle(pca_file)
        input_data_index = list(X.index)

    else:
        print("read from raw file")
        input_data = sc.read_h5ad(input_file)
        input_data_index = list(input_data.obs.index)
        X = pd.DataFrame(input_data.obsm[obsm_key].copy(),index=input_data_index)
        
        X.to_pickle(pca_file)
        del input_data


    print("read pickle")
    if os.path.isfile(dict_file):
        print("read from dictionary")
        with open(dict_file, mode='br') as fi:
            gRNA_dict = pickle.load(fi)
    else:
        print("read from raw count")
        sgRNA_data = pd.read_pickle(sgRNA_file)
        sgRNA_data = sgRNA_data.T
        sgRNA_data = sgRNA_data[sgRNA_data.index.isin(input_data_index)]
        
        print("finish reading raw count")
        non_zero_location = np.where(sgRNA_data.values != 0)
        gRNA_name_list = sgRNA_data.columns.to_list()
        cell_name_list = sgRNA_data.index.to_list()

        gRNA_convert = [gRNA_name_list[i] for i in non_zero_location[1]]
        cell_convert = [cell_name_list[i] for i in non_zero_location[0]]
        non_zero_location_zip = list(zip(gRNA_convert,cell_convert))
        non_zero_array = np.array(non_zero_location_zip)
        
        gRNA_idx_list = np.unique(non_zero_array[:,0])
        gRNA_dict = {}
        for gRNA_name in tqdm(gRNA_idx_list):
            gRNA_dict[gRNA_name] = non_zero_array[non_zero_array[:,0]==gRNA_name][:,1].tolist()
        
        with open(dict_file, mode='wb') as fo:
            pickle.dump(gRNA_dict, fo)
        
        del sgRNA_data
        gc.collect()
    gc.collect()
    return (X,gRNA_dict)

def pairwise_torch(X,cell_id_list,device,batch_size=100,norm_p=2.0):
    cell_num_len = np.array([len(x) for x in cell_id_list])
    max_cell_num = max(cell_num_len)
    res = []

    test_tensor_nonpad = []
    test_tensor = []
    musk_tensor = []
    for index,cell_names in enumerate(cell_id_list):
        pad_num = max_cell_num - cell_num_len[index]
        test_tensor += [torch.tensor(np.pad(X.loc[cell_names,:],((pad_num,0),(0,0))))]
        musk_tensor += [torch.tensor([[0]*pad_num+[1]*cell_num_len[index]])]
    test_tensor = torch.stack(test_tensor).to(device)
    musk_tensor = torch.stack(musk_tensor).to(device)

    res = []
    cell_num_len_tensor = torch.tensor(cell_num_len).to(device)
    
    combis = np.array(list(combinations(range(len(cell_id_list)), 2)) 
                      + [(x,x) for x in range(len(cell_id_list))])
    total_len_combis = np.arange(len(combis))
    calc_list = np.array_split(total_len_combis,len(combis) // batch_size)

    for combis_idx in tqdm(calc_list):
        batch_sample_size = len(combis_idx)
        tensor_0 = test_tensor[combis[combis_idx][:,0]]
        tensor_1 = test_tensor[combis[combis_idx][:,1]]

        musk_2D = torch.mul(
            musk_tensor[combis[combis_idx][:,0]].view(batch_sample_size,-1,1),
            musk_tensor[combis[combis_idx][:,1]]
        )
        sum_cross = torch.cdist(tensor_0,tensor_1,p=norm_p).pow(2).mul(musk_2D).sum((1,2))
        sum_cross = sum_cross.mul((cell_num_len_tensor[combis[combis_idx][:,0]] 
                                   * cell_num_len_tensor[combis[combis_idx][:,1]]).reciprocal())
        res += [sum_cross]
    combi_1 = combis[:,0]
    combi_2 = combis[:,1]
    
    return_res = torch.concat(res).to("cpu").tolist()
    
    del test_tensor,musk_tensor,cell_num_len_tensor
    gc.collect()
    torch.cuda.empty_cache
    
    return list(zip(combi_1,combi_2,return_res))


def permutation_test(X,test_cell1,test_cell2,device,batch_num=10,total_permute=1000,
                     norm_p=2.0,return_permute=True):
    test_cells_concat = np.concatenate((test_cell1, test_cell2)).tolist()
    repeat_num = int(total_permute/batch_num)
    num_cell1 = len(test_cell1)
    num_cell2 = len(test_cell2)
    num_all = num_cell1 + num_cell2
    
    all_cells_dist = torch.Tensor(X.loc[test_cells_concat,:].to_numpy()).to(device)
    pairwise_dist_tmp = []
    
    for index_arr in np.array_split(np.arange(num_all),1):
        pairwise_dist_tmp += [torch.cdist(all_cells_dist[index_arr],all_cells_dist,p=norm_p).pow(2)]

    pairwise_dist = torch.cat(pairwise_dist_tmp).repeat(batch_num,1,1)
    
    sum_target = pairwise_dist[0][num_cell1:][:,num_cell1:].sum() / (num_cell2*num_cell2)
    sum_non_target = pairwise_dist[0][:num_cell1][:,:num_cell1].sum() / (num_cell1*num_cell1)
    sum_cross = pairwise_dist[0][:num_cell1][:,num_cell1:].sum() / (num_cell1*num_cell2)

    obs_edist = 2*sum_cross-sum_non_target-sum_target
    if return_permute==False:
        return obs_edist
    else:
        e_dist_list = []

        np.random.seed(0)

        for i in range(repeat_num):
            random_id = np.array([np.random.permutation(num_all) for i in range(batch_num)])

            group_target = torch.tensor(random_id[:,:num_cell1]).to(device)
            group_non_target = torch.tensor(random_id[:,num_cell1:]).to(device)

            extracted = pairwise_dist.gather(1,group_target.unsqueeze(2).expand(-1, -1, num_all))
            extracted = extracted.gather(2, group_target.unsqueeze(1).expand(-1, num_cell1, -1))
            sum_target = extracted.sum((1,2)) / (num_cell1*num_cell1)

            extracted = pairwise_dist.gather(1,group_non_target.unsqueeze(2).expand(-1, -1, num_all))
            extracted = extracted.gather(2, group_non_target.unsqueeze(1).expand(-1, num_cell2, -1))
            sum_non_target = extracted.sum((1,2)) / (num_cell2*num_cell2)

            extracted = pairwise_dist.gather(1,group_target.unsqueeze(2).expand(-1, -1, num_all))
            extracted = extracted.gather(2, group_non_target.unsqueeze(1).expand(-1, num_cell1, -1))
            sum_cross = extracted.sum((1,2)) / (num_cell2*num_cell1)

            e_dist_list += [2*sum_cross-sum_target-sum_non_target]

        e_dist_list = torch.cat(e_dist_list)
        
        #Move to cpu memory
        e_dist_list = e_dist_list.to("cpu")
        obs_edist = obs_edist.to("cpu")
        
        #cleaning to release GPU memory
        del pairwise_dist_tmp,sum_target,sum_non_target,sum_cross
        gc.collect()
        torch.cuda.empty_cache

        return (obs_edist,e_dist_list)


def disco_test(X,test_cell_list,device,batch_num=5,total_permute=1000,norm_p=2.0):
    alpha=1
    test_cells_concat = np.concatenate(test_cell_list).tolist()
    num_of_samples = len(test_cell_list)
    
    length_sample_array = [len(i) for i in test_cell_list]
    sample_index_array = [[i]*length_sample_array[i] for i in range(num_of_samples)]
    sample_index_array = np.concatenate(sample_index_array)
    num_all = len(test_cells_concat)
    
    repeat_num = int(total_permute/batch_num)
    all_cells_dist = torch.Tensor(X.loc[test_cells_concat,:].to_numpy()).to(device)
    pairwise_dist_tmp = []
    
    for index_arr in np.array_split(np.arange(num_all),1):
        pairwise_dist_tmp += [torch.cdist(all_cells_dist[index_arr],all_cells_dist,p=norm_p).pow(alpha)]

    pairwise_dist = torch.cat(pairwise_dist_tmp).repeat(batch_num,1,1)
    total_disp = (num_all/2)*(pairwise_dist[0].sum()/num_all/num_all)
    
    within_disp_list = []
    for i in range(num_of_samples):
        cell_index = (sample_index_array == i)
        tmp = (length_sample_array[i]/2)*pairwise_dist[0][cell_index,:][:,cell_index].sum()/length_sample_array[i]/length_sample_array[i]
        within_disp_list.append(tmp)
    within_disp = torch.stack(within_disp_list).sum()
    
    between_disp = total_disp - within_disp
    F_value = (between_disp/(num_of_samples-1))/(within_disp/(num_all-num_of_samples))

    np.random.seed(0)
    F_value_permute_list = []
    
    for i in range(repeat_num):
        random_id = np.array([np.random.permutation(num_all) for i in range(batch_num)])
        
        group_cell_id = []
        for k in range(num_of_samples):  
            cell_index = (sample_index_array == k)
            group_cell_id.append(torch.tensor(random_id[:,cell_index]).to(device))
        
        within_disp_list = []
        for k in range(num_of_samples):
            extracted = pairwise_dist.gather(1,group_cell_id[k].unsqueeze(2).expand(-1, -1, num_all))
            extracted = extracted.gather(2, group_cell_id[k].unsqueeze(1).expand(-1, length_sample_array[k], -1))
            extracted = (length_sample_array[k]/2)*extracted
            within_disp_list.append(extracted.sum((1,2)) / (length_sample_array[k]*length_sample_array[k]))
        within_disp_permute = torch.sum(torch.stack(within_disp_list),0)
        total_disp_permute = total_disp.expand(batch_num)
        between_disp_permute = total_disp_permute - within_disp_permute
        
        F_value_permute_list += [(between_disp_permute/(num_of_samples-1))/(within_disp_permute/(num_all-num_of_samples))]
    F_value_permute_list = torch.cat(F_value_permute_list).to("cpu")
    
    del all_cells_dist,between_disp_permute,within_disp_permute,total_disp_permute,pairwise_dist
    gc.collect()
    torch.cuda.empty_cache
        
    return (F_value,F_value_permute_list)

def get_unique_list(seq):
    return_list = []
    seen = []
    for x in seq:
        if x in seen:
            continue
        elif (x[1],x[0]) in seen:
            continue
        else:
            seen.append(x)
    return seen

def extract_gene_name(entry):
    if entry.startswith('OR'):
        return entry.split("-")[0]
    else:
        return entry.split("_")[0]
    
def extract_transcript_name(entry):
    gene_name = ""
    transcript_name = ""
    if entry.startswith('OR'):
        gene_name = entry.split("-")[0]
    else:
        gene_name = entry.split('_')[0]
        transcript_name = entry.split('-')[-2:][0]
    if len(entry.split('_')) <3:
        transcript_name = gene_name
    else:
        transcript_name = gene_name + ':' + transcript_name
    return transcript_name