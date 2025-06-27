# calculate dssp for each PPI structural model and save as separate files
import multiprocessing
import os
import pickle
import pandas as pd
import subprocess
from collections import defaultdict
from lib_sasa_scan import load_structure, sasa_scan
from pathlib import Path
from typing import Union

def default_str_dict():
    return defaultdict(str)

def default_list_dict():
    return defaultdict(list)

aaInfo_all = defaultdict(default_str_dict)
RSAInfo_all = defaultdict(default_list_dict)
deltaRSAInfo_all = defaultdict(default_list_dict)

def round_delta_RSA(residue: dict[str, float], 
                    threshold: float) -> dict[str, Union[str, float]]:
    '''
    Set delta RSA to 0 if it is below a specified threshold. 
    '''
    if residue['dsasa'] >= threshold or residue['dsasa'] == 0.:
        return residue
    else:
        return {'aa': residue['aa'], 'rsaMonomer': residue['rsaMonomer'],'rsaComplex':residue['rsaComplex'] ,'dsasa': 0.0}

def calculate_rsa_protein(complexID: str, 
                          proteinID: str, 
                          chainID: str, 
                          partnerChain: str, 
                          ppiDir: Union[str, Path], 
                          path_to_dssp: Union[str, Path], 
                          dataset='strint') -> None: 
    '''
    Calculate RSA for a protein. 

    Args:
        complexID (str): Complex ID where PPI is. 
        proteinID (str): Protein name. 
        chainID (str): Chain ID of the protein in the PPI. 
        partnerChain (str): Chain ID of the other protein in the PPI. 
        ppiDir (Union[str, Path]): Directory of PPI structural models. 
        path_to_dssp (Union[str, Path]): Path to executable dssp. 
        dataset (str): Either strint or alphafold-resettafold. 

    Returns:
        None
    '''
    global aaInfo_all, RSAInfo_all, deltaRSAInfo_all
    minNonZeroDeltaRSA = 0.001
    if dataset == 'strint': 
        struct = load_structure(ppiDir / ('pdb' + complexID + '.ent'))
    elif dataset == 'alphafold-rosettafold':
        p1, p2 = complexID.split('=')
        complexID_af = '_'.join([p1, p2])
        struct = load_structure(ppiDir / (complexID_af + '.pdb'))
    else: 
        print("Please set 'dataset' to either 'strint' or 'alphafold-resettafold'. ")
        return
    sasa = sasa_scan(struct, chainID, chainID + partnerChain, path_to_dssp, True, tmpPath='./atom.pdb')[0]
    resInfo = dict([(k, round_delta_RSA(v, minNonZeroDeltaRSA)) for k, v in sasa.items()])
    for pos, info_d in resInfo.items():
        aaInfo_all[proteinID][pos] = info_d['aa']
        RSAInfo_all[proteinID][pos].append(info_d['sasa'])
        deltaRSAInfo_all[proteinID][pos].append(info_d['dsasa'])

def calculate_rsa_complex(complexID, ppiDir, path_to_dssp, dataset='strint'): 
    '''
    Calculate RSA for each protein in a PPI. 

    Args:
        complexID (str): Complex ID where PPI is. 
        ppiDir (Union[str, Path]): Directory of PPI structural models. 
        path_to_dssp (Union[str, Path]): Path to executable dssp. 
        dataset (str): Either strint or alphafold-resettafold. 
    
    Returns:
        None
    '''
    global aaInfo_all, RSAInfo_all, deltaRSAInfo_all
    print(complexID)
    p1, p2 = complexID.split('=')
    calculate_rsa_protein(complexID, p1, 'A', 'B', ppiDir, path_to_dssp, dataset=dataset)
    calculate_rsa_protein(complexID, p2, 'B', 'A', ppiDir, path_to_dssp, dataset=dataset)

def main():

    global aaInfo_all, RSAInfo_all, deltaRSAInfo_all

    dataDir = Path('../data')

    extDir = dataDir / 'external'

    procDir = dataDir / 'processed'

    afDir = procDir / 'AlphaFold2-RoseTTAFold'

    homoDir = procDir / 'ppi_models'

    # input files

    strIntFile = procDir / 'structural_interactome_merged.txt'

    # output files
    RSAFile = procDir / 'RSA_dicts.pkl'

    path_to_dssp = 'mkdssp'

    strInt = pd.read_table(strIntFile)

    for index, row in strInt.iterrows():
        if row['Template_file_ID'] == 'AlphaFold2':
            complexID = row['Chain_pairs'].split('_')[0]
            fileName = f'{complexID}.pdb'
            ppiDir = afDir
            dataset = 'alphafold-rosettafold'
            calculate_rsa_complex(complexID, ppiDir, path_to_dssp, dataset=dataset)
        else:
            complexID = row['Complex_ID']
            fileName = f'pdb{complexID}.ent'
            ppiDir = homoDir
            dataset = 'strint'
            calculate_rsa_complex(complexID, ppiDir, path_to_dssp, dataset=dataset)
    pickle.dump([aaInfo_all, RSAInfo_all, deltaRSAInfo_all], open(RSAFile, 'wb'))

if __name__ == '__main__':
    main()