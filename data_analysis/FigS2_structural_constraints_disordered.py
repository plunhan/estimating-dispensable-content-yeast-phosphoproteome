'''
Categorize conditional p-sites, universal p-sites, and random S/T according to the exposure
and compare the evolutionary rate in each exposure category. 
'''

import pandas as pd
import pickle
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path
from plot_tools import plot_consurf_exposure
from proteomic_tools import (retrieve_references_by_residue_type, 
							 get_phosphosites_given_perturbations,
							 sample_random_sites, 
							 retrieve_references_by_order, 
                             retrieve_exposure_references, 
							 calculate_exposure_consurf_alphafold, 
                             calculate_exposed_relative)
from scipy.stats import ranksums

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
    #sample_residues = 'ACDEFGHIKLMNPQRSTVWY'

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
    rsa_pkl = procDir / 'RSA_dicts.pkl'
    lanz90PKL = paperDir / 'lanz90.pkl'
    rsa_alphafold_pkl = paperDir / 'rsa_alphafold.pkl'

    # output files
    Fig6A = paperDir / f'Figure 6A.jpg'

    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    ultradeep = pickle.load(open(ultradeepPKL, 'rb'))
    phosStres = pickle.load(open(phosStresPKL, 'rb'))
    consurf = pickle.load(open(consurfPKL, 'rb'))
    diso = pickle.load(open(disoPKL, 'rb'))
    sequences = pickle.load(open(sequencePKL, 'rb'))
    sgd = pickle.load(open(sgdPKL, 'rb'))
    biogrid = pickle.load(open(biogridPKL, 'rb'))
    lanz90 = pickle.load(open(lanz90PKL, 'rb'))
    aaInfo, RSAInfo, dRSAInfo = parse_rsa_d(rsa_pkl, IDMappingDict, method=method)
    rsa_alphafold = pickle.load(open(rsa_alphafold_pkl, 'rb'))

    phosStres = {key: retrieve_references_by_residue_type(references, sample_residues) for key, references in phosStres.items()}
    cond_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 11)))
    univ_psites = get_phosphosites_given_perturbations(phosStres, list(range(92, 102)))
    all_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 102)))
    cond_psites_ST = retrieve_references_by_residue_type(cond_psites, sample_residues)
    univ_psites_ST = retrieve_references_by_residue_type(univ_psites, sample_residues)
    all_psites_ST = retrieve_references_by_residue_type(all_psites, sample_residues)
    cond_psites_dis = retrieve_references_by_order(cond_psites_ST, diso, 'disordered')
    univ_psites_dis = retrieve_references_by_order(univ_psites_ST, diso, 'disordered')
    all_psites_dis = retrieve_references_by_order(all_psites_ST, diso, 'disordered')
    cond_psites_ord = retrieve_references_by_order(cond_psites_ST, diso, 'ordered')
    univ_psites_ord = retrieve_references_by_order(univ_psites_ST, diso, 'ordered')
    all_psites_ord = retrieve_references_by_order(all_psites_ST, diso, 'ordered')

    exclusions = ultradeep.union(sgd, biogrid, lanz90) # All reported p-sites
    randomST = sample_random_sites(ultradeep, exclusions, sequences, sample_residues)
    randomST_dis = retrieve_references_by_order(randomST, diso, 'disordered')
    randomST_ord = retrieve_references_by_order(randomST, diso, 'ordered')

    # Calculate p-values
    randomST_dis_exposed = calculate_exposed_relative(randomST_dis, consurf, rsa_alphafold, dRSAInfo)
    univ_dis_exposed = calculate_exposed_relative(univ_psites_dis, consurf, rsa_alphafold, dRSAInfo)
    cond_dis_exposed = calculate_exposed_relative(cond_psites_dis, consurf, rsa_alphafold, dRSAInfo)

    _, pvalue = ranksums(univ_dis_exposed, cond_dis_exposed)
    print(f'univ vs cond: {pvalue}')
    _, pvalue = ranksums(randomST_dis_exposed, cond_dis_exposed)
    print(f'rand vs cond: {pvalue}')
    _, pvalue = ranksums(univ_dis_exposed, randomST_dis_exposed)
    print(f'univ vs rand: {pvalue}')

    randomST_ord_exposed = calculate_exposed_relative(randomST_ord, consurf, rsa_alphafold, dRSAInfo)
    univ_ord_exposed = calculate_exposed_relative(univ_psites_ord, consurf, rsa_alphafold, dRSAInfo)
    cond_ord_exposed = calculate_exposed_relative(cond_psites_ord, consurf, rsa_alphafold, dRSAInfo)

    _, pvalue = ranksums(univ_ord_exposed, cond_ord_exposed)
    print(f'univ vs cond: {pvalue}')
    _, pvalue = ranksums(randomST_ord_exposed, cond_ord_exposed)
    print(f'rand vs cond: {pvalue}')
    _, pvalue = ranksums(univ_ord_exposed, randomST_ord_exposed)
    print(f'univ vs rand: {pvalue}')

    # Plot figure
    df_ls = []
    for i, references in enumerate([randomST_dis, cond_psites_dis, univ_psites_dis]): 
        if i == 0: 
            resType = "Non-phosphorylated S/T sites"
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
    order = ['Non-phosphorylated S/T sites', 'Conditional phosphosites', 'Universal phosphosites']
    plot_consurf_exposure(df, order, Fig6A, figFmt, 'disordered', (-0.5, 0.3))

if __name__ == '__main__':
    main()