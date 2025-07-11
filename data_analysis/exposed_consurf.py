'''
Categorize conditional p-sites, universal p-sites, and random S/T according to the exposure
and compare the evolutionary rate in each exposure category. 
'''

import numpy as np
import pandas as pd
import pickle
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path
from plot_tools import plot_consurf_exposure
from proteomic_tools import (retrieve_references_by_residue_type, 
							 get_phosphosites_given_perturbations,
                             retrieve_ConSurf_score, 
							 sample_random_sites, 
							 retrieve_references_by_order, 
							 calculate_exposure_consurf)

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

def calculate_exposed_consurf(exposed_residues_consurf_d: dict[str, float], 
                              sites_d: dict[str, set[str]],
                              consurf_d: dict[str, dict[int, float]], 
                              diso: dict[str, dict[int, float]],
                              aaInfo: dict[str, dict[int, str]], 
                              RSAInfo: dict[str, dict[int, float]], 
                              dRSAInfo: dict[str, dict[int, float]]) -> dict[str, float]:
    for key, sites in sites_d.items(): 
        if key == 'randomAll':
            sample_residues = 'ACDEFGHIKLMNPQRSTVWY'
        else:
            sample_residues = 'ST'
        interfacial_psite = []
        exposed_psite = []
        buried_psite = []
        psites_ST = retrieve_references_by_residue_type(sites, sample_residues)
        psites_dis = retrieve_references_by_order(psites_ST, diso, 'disordered')
        for psite in psites_dis: 
            systematic_name, site = psite.split('_')
            if systematic_name in ['YKL021C', 'YDR098C', 'YMR112C']:
                continue
            aa = site[0]
            position = int(site[1:])
            try: 
                if dRSAInfo[systematic_name][position] > 0:
                    interfacial_psite.append(psite)
                elif RSAInfo[systematic_name][position] > 0.25:
                    exposed_psite.append(psite)
                elif RSAInfo[systematic_name][position] <= 0.25:
                    buried_psite.append(psite)
            except KeyError:
                continue
        _, consurf_exposed = retrieve_ConSurf_score(exposed_psite, consurf_d)
        exposed_residues_consurf_d[key] = np.median(consurf_exposed)
    return exposed_residues_consurf_d

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
    rsa_pkl = procDir / 'RSA_dicts.pkl'

    # output files
    Fig6 = paperDir / f'Figure 6 {method}.jpg'

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

    phosStres = {key: retrieve_references_by_residue_type(references, sample_residues) for key, references in phosStres.items()}
    cond_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 11)))
    univ_psites = get_phosphosites_given_perturbations(phosStres, list(range(92, 102)))
    all_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 102)))

    exposed_consurf_d = {}

    exclusions = ultradeep.union(sgd, biogrid) # All reported p-sites
    randomST = sample_random_sites(ultradeep, ultradeep, sequences, sample_residues)
    randomST_all = sample_random_sites(ultradeep, exclusions, sequences, sample_residues)
    randomAll = sample_random_sites(ultradeep, exclusions, sequences, 'ACDEFGHIKLMNPQRSTVWY')

    sites_d = {
        'cond': cond_psites, 
        'univ': univ_psites, 
        'all': all_psites, 
        'biogrid': biogrid,
        'sgd': sgd, 
        'lanz90': lanz90, 
        'randomST': randomST,
        'randomST_all': randomST_all,
        'randomAll': randomAll
    }

    exposed_consurf_d = calculate_exposed_consurf(exposed_consurf_d, 
                                                  sites_d,
                                                  consurf, 
                                                  diso, 
                                                  aaInfo, 
                                                  RSAInfo, 
                                                  dRSAInfo)

    for key, consurf_median in exposed_consurf_d.items():
        print(f'{key}:\t{consurf_median}')

if __name__ == '__main__':
    main()