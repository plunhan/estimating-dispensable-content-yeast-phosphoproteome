'''
Estimate fraction of functionally dispensable p-sites in phosphoproteomes. 
'''

import numpy as np
import os
import pandas as pd
import pickle
from estimate_tools import estimate_pi_mixture_model
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path
from plot_tools import plot_consurf_exposure
from proteomic_tools import (get_phosphosites_given_perturbations, 
							 retrieve_references_by_order, 
							 retrieve_references_by_residue_type, 
							 retrieve_ConSurf_score,
                             retrieve_exposure_references, 
							 sample_random_sites, 
                             calculate_exposure_consurf_alphafold)
from scipy.stats import ranksums

def retrieve_exposure_references_alphafold(references, rsa_alphafold, dRSAInfo):
    result_set = set()
    for reference in references: 
        systematic_name, site = reference.split('_')
        if systematic_name in ['YKL021C', 'YDR098C', 'YMR112C']:
            continue
        aa = site[0]
        position = int(site[1:])
        try: 
            if rsa_alphafold[systematic_name][position] > 0.25 and dRSAInfo[systematic_name][position] == 0:
                result_set.add(reference)
        except KeyError:
            continue
    return result_set

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
    rsa_alphafold_pkl = paperDir / 'rsa_alphafold.pkl'

    # output files
    FigS2C = paperDir / 'Figure S2C.jpg'

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
    rsa_alphafold = pickle.load(open(rsa_alphafold_pkl, 'rb'))

    # Part I: difference between conditional and universal phosphosites
    cond_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 11)))
    univ_psites = get_phosphosites_given_perturbations(phosStres, list(range(92, 102)))
    all_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 102)))

    cond_psites_ST = retrieve_references_by_residue_type(cond_psites, sample_residues)
    univ_psites_ST = retrieve_references_by_residue_type(univ_psites, sample_residues)
    all_psites_ST = retrieve_references_by_residue_type(all_psites, sample_residues)

    cond_psites_ord = retrieve_references_by_order(cond_psites_ST, diso, 'ordered')
    univ_psites_ord = retrieve_references_by_order(univ_psites_ST, diso, 'ordered')
    all_psites_ord = retrieve_references_by_order(all_psites_ST, diso, 'ordered')

    cond_psites_ord_exp = retrieve_exposure_references_alphafold(cond_psites_ord, rsa_alphafold, dRSAInfo)
    univ_psites_ord_exp = retrieve_exposure_references_alphafold(univ_psites_ord, rsa_alphafold, dRSAInfo)

    cond_ord_consurf_references, cond_ord_consurf = retrieve_ConSurf_score(cond_psites_ord_exp, consurf)
    univ_ord_consurf_references, univ_ord_consurf = retrieve_ConSurf_score(univ_psites_ord_exp, consurf)

    exclusions = ultradeep.union(sgd, biogrid, lanz90) # All reported p-sites
    randomST = sample_random_sites(ultradeep, exclusions, sequences, sample_residues)
    randomST_ord = retrieve_references_by_order(randomST, diso, 'ordered')
    randomST_psites_ord_exp = retrieve_exposure_references_alphafold(randomST_ord, rsa_alphafold, dRSAInfo)
    randomST_ord_consurf_references, randomST_ord_consurf = retrieve_ConSurf_score(randomST_psites_ord_exp, consurf)

    disordered_data = [cond_ord_consurf, randomST_ord_consurf, univ_ord_consurf]
    labels = ['Non-phosphorylated S/T', 'Conditional', 'Universal']

    print(np.median(cond_ord_consurf))
    print(np.median(univ_ord_consurf))
    print(np.median(randomST_ord_consurf))
    print(ranksums(cond_ord_consurf, univ_ord_consurf))
    print(ranksums(randomST_ord_consurf, univ_ord_consurf))
    print(ranksums(cond_ord_consurf, randomST_ord_consurf))

    df_ls = []
    for i, references in enumerate([randomST_ord, cond_psites_ord, univ_psites_ord]): 
        if i == 0: 
            resType = "Non-phosphorylated S/T"
        elif i == 1: 
            resType = "Conditional phosphosites"
        else: 
            resType = "Universal phosphosites"
        rows = calculate_exposure_consurf_alphafold(references, resType, consurf, rsa_alphafold, dRSAInfo)
        if i == 1:
            df_ls = rows + df_ls
        else:
            df_ls.extend(rows)

    df = pd.DataFrame(df_ls, columns=['Exposure', 'Type', 'Median', 'Standard error'])
    print(df)
    order = ['Non-phosphorylated S/T', 'Conditional phosphosites', 'Universal phosphosites']
    plot_consurf_exposure(df, order, FigS2C, figFmt, (-0.2, 0.2))

if __name__ == '__main__':
    main()