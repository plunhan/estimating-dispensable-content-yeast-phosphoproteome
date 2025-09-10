'''
Estimate fraction of functionally dispensable p-sites in phosphoproteomes. 
'''

import numpy as np
import os
import pickle
from estimate_tools import estimate_pi_mixture_model
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path
from proteomic_tools import (get_phosphosites_given_perturbations, 
							 retrieve_references_by_order, 
							 retrieve_references_by_residue_type, 
							 retrieve_ConSurf_score,
                             retrieve_exposure_references, 
							 sample_random_sites)

def default_list_dict():
    return defaultdict(list)

def default_str_dict():
    return defaultdict(str)

def parse_rsa_d(inPath, IDMappingDict, method='maximum'):
    # method can be average or maximum
    aaInfo_all, RSAInfo_all, deltaRSAInfo_all = pickle.load(open(inPath, 'rb'))
    RSAInfo_filtered = {}
    if method == "maximum":
        for protein, residues in RSAInfo_all.items():
            RSAInfo_filtered[protein] = {resnum: max(RSAs) for resnum, RSAs in residues.items()}
    elif method == "average":
        for protein, residues in RSAInfo_all.items():
            RSAInfo_filtered[protein] = {resnum: sum(RSAs) / len(RSAs) for resnum, RSAs in residues.items()}
    elif method == "minimum":
        for protein, residues in RSAInfo_all.items():
            RSAInfo_filtered[protein] = {resnum: min(RSAs) for resnum, RSAs in residues.items()}
    deltaRSAInfo_filtered = {}
    for protein, residues in deltaRSAInfo_all.items():
        deltaRSAInfo_filtered[protein] = {resnum: any(deltaRSA > 0 for deltaRSA in deltaRSAs) for resnum, deltaRSAs in residues.items()}
    aaInfo_all = {IDMappingDict[key]: value for key, value in aaInfo_all.items() if key in IDMappingDict}
    RSAInfo_filtered = {IDMappingDict[key]: value for key, value in RSAInfo_filtered.items() if key in IDMappingDict}
    deltaRSAInfo_filtered = {IDMappingDict[key]: value for key, value in deltaRSAInfo_filtered.items() if key in IDMappingDict}
    return aaInfo_all, RSAInfo_filtered, deltaRSAInfo_filtered

def main():

    method = 'maximum'

    figFmt = 'jpg'

    sample_residues = 'ST'

    dataDir = Path('../../data')

    extDir = dataDir / 'external'

    procDir = dataDir / 'processed'

    paperDir = procDir / 'paper'

    # input files
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    ultradeepPKL = paperDir / 'ultradeep_reference_phosphoproteome.pkl'
    phosStresPKL = paperDir / 'quantitative_phosphosites.pkl'
    consurfPKL = paperDir / 'consurf_all.pkl'
    disoPKL = paperDir / 'diso_all.pkl'
    sequencePKL = paperDir / 'Scer_seq.pkl'
    sgdPKL = paperDir / 'SGD.pkl'
    biogridPKL = paperDir / 'BioGRID.pkl'
    lanz90PKL = paperDir / 'lanz90.pkl'
    lanz70PKL = paperDir / 'lanz70.pkl'
    rsa_pkl = procDir / 'RSA_dicts.pkl'

    # output files
    Fig2A = paperDir / 'Figure 2A.jpg'
    Fig2B = paperDir / 'Figure 2B.jpg'

    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    ultradeep = pickle.load(open(ultradeepPKL, 'rb'))
    phosStres = pickle.load(open(phosStresPKL, 'rb'))
    consurf = pickle.load(open(consurfPKL, 'rb'))
    diso = pickle.load(open(disoPKL, 'rb'))
    sequences = pickle.load(open(sequencePKL, 'rb'))
    sgd = pickle.load(open(sgdPKL, 'rb'))
    biogrid = pickle.load(open(biogridPKL, 'rb'))
    lanz90 = pickle.load(open(lanz90PKL, 'rb'))
    lanz70 = pickle.load(open(lanz70PKL, 'rb'))
    aaInfo, RSAInfo, dRSAInfo = parse_rsa_d(rsa_pkl, IDMappingDict, method=method)

    # Part I: difference between conditional and universal phosphosites
    cond_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 11)))
    univ_psites = get_phosphosites_given_perturbations(phosStres, list(range(92, 102)))
    all_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 102)))

    cond_psites_ST = retrieve_references_by_residue_type(cond_psites, sample_residues)
    univ_psites_ST = retrieve_references_by_residue_type(univ_psites, sample_residues)
    all_psites_ST = retrieve_references_by_residue_type(all_psites, sample_residues)

    cond_psites_dis = retrieve_references_by_order(cond_psites_ST, diso, 'disordered')
    univ_psites_dis = retrieve_references_by_order(univ_psites_ST, diso, 'disordered')
    all_psites_dis = retrieve_references_by_order(all_psites_ST, diso, 'disordered')

    cond_psites_dis_exp = retrieve_exposure_references(cond_psites_dis, sample_residues, aaInfo, RSAInfo, dRSAInfo)
    univ_psites_dis_exp = retrieve_exposure_references(univ_psites_dis, sample_residues, aaInfo, RSAInfo, dRSAInfo)

    cond_dis_consurf_references, cond_dis_consurf = retrieve_ConSurf_score(cond_psites_dis_exp, consurf)
    univ_dis_consurf_references, univ_dis_consurf = retrieve_ConSurf_score(univ_psites_dis_exp, consurf)

    # exclusions = ultradeep
    exclusions = ultradeep.union(sgd, biogrid, lanz90) # All reported p-sites
    # exclusions = set()
    randomST = sample_random_sites(ultradeep, exclusions, sequences, sample_residues)
    randomST_dis = retrieve_references_by_order(randomST, diso, 'disordered')
    randomST_psites_dis_exp = retrieve_exposure_references(randomST_dis, sample_residues, aaInfo, RSAInfo, dRSAInfo)
    randomST_dis_consurf_references, randomST_dis_consurf = retrieve_ConSurf_score(randomST_psites_dis_exp, consurf)

    disordered_data = [cond_dis_consurf, randomST_dis_consurf, univ_dis_consurf]
    labels = ['Random S/T', 'Conditional', 'Universal']

    print(np.median(cond_dis_consurf))
    print(np.median(univ_dis_consurf))
    print(np.median(randomST_dis_consurf))
    print(len(randomST_dis_consurf_references), len(randomST_dis_consurf))

    # for db in ['sgd', 'biogrid', 'lanz70', 'lanz90', 'leutert']:
    #     if db == 'sgd': 
    #         psites = sgd
    #     elif db == 'biogrid':
    #         psites = biogrid
    #     elif db == 'lanz70':
    #         psites = lanz70
    #     elif db == 'lanz90':
    #         psites = lanz90
    #     elif db == 'leutert':
    #         psites = all_psites
    #     print(db)
    #     all_psites_ST_db = retrieve_references_by_residue_type(psites, sample_residues)
    #     all_psites_dis_db = retrieve_references_by_order(all_psites_ST_db, diso, 'disordered')
    #     all_psites_dis_exp = retrieve_exposure_references(all_psites_dis_db, sample_residues, aaInfo, RSAInfo, dRSAInfo)
    #     all_dis_consurf_references_db, all_dis_consurf_db = retrieve_ConSurf_score(all_psites_dis_exp, consurf)
    #     print(f'{np.median(all_dis_consurf_db):.3f} of {len(all_dis_consurf_db)}')
    #     dispensable = estimate_pi_mixture_model(cond_dis_consurf, univ_dis_consurf, all_dis_consurf_db)
    #     print(f'{db}, {dispensable:.2f}, {(np.median(all_dis_consurf_db) - np.median(univ_dis_consurf))/(np.median(cond_dis_consurf) - np.median(univ_dis_consurf)):.2f}')

if __name__ == '__main__':
    main()