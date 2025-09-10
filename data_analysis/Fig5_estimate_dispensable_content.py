'''
Estimate fraction of functionally dispensable p-sites in phosphoproteomes. 
'''

import numpy as np
import os
import pickle
from estimate_tools import estimate_pi_mixture_model
from pathlib import Path
from proteomic_tools import (get_phosphosites_given_perturbations, 
							 retrieve_references_by_order, 
							 retrieve_references_by_residue_type, 
							 retrieve_ConSurf_score,
							 sample_random_sites, 
                             calculate_consurf_difference_psites)
from scipy.stats import ranksums

def main():

    figFmt = 'jpg'

    sample_residues = 'ST'

    dataDir = Path('../../data')

    extDir = dataDir / 'external'

    procDir = dataDir / 'processed'

    paperDir = procDir / 'paper'

    # input files
    ultradeepPKL = paperDir / 'ultradeep_reference_phosphoproteome.pkl'
    phosStresPKL = paperDir / 'quantitative_phosphosites.pkl'
    consurfPKL = paperDir / 'consurf_all.pkl'
    disoPKL = paperDir / 'diso_all.pkl'
    sequencePKL = paperDir / 'Scer_seq.pkl'
    sgdPKL = paperDir / 'SGD.pkl'
    biogridPKL = paperDir / 'BioGRID.pkl'
    lanz90PKL = paperDir / 'lanz90.pkl'
    lanz70PKL = paperDir / 'lanz70.pkl'
    PhosphoGIRD_PKL = paperDir / 'PhosphoGIRD.pkl'

    # output files
    Fig2A = paperDir / 'Figure 2A.jpg'
    Fig2B = paperDir / 'Figure 2B.jpg'

    ultradeep = pickle.load(open(ultradeepPKL, 'rb'))
    phosStres = pickle.load(open(phosStresPKL, 'rb'))
    consurf = pickle.load(open(consurfPKL, 'rb'))
    diso = pickle.load(open(disoPKL, 'rb'))
    sequences = pickle.load(open(sequencePKL, 'rb'))
    sgd = pickle.load(open(sgdPKL, 'rb'))
    biogrid = pickle.load(open(biogridPKL, 'rb'))
    lanz90 = pickle.load(open(lanz90PKL, 'rb'))
    lanz70 = pickle.load(open(lanz70PKL, 'rb'))
    PhosphoGIRD = pickle.load(open(PhosphoGIRD_PKL, 'rb'))

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

    cond_relative = calculate_consurf_difference_psites(cond_psites_dis, consurf, window_size=5)
    univ_relative = calculate_consurf_difference_psites(univ_psites_dis, consurf, window_size=5)

    cond_dis_consurf_references, cond_dis_consurf = retrieve_ConSurf_score(cond_psites_dis, consurf)
    univ_dis_consurf_references, univ_dis_consurf = retrieve_ConSurf_score(univ_psites_dis, consurf)

    exclusions = ultradeep.union(sgd, biogrid, lanz90) # All reported p-sites
    randomST = sample_random_sites(ultradeep, exclusions, sequences, sample_residues)
    randomST_dis = retrieve_references_by_order(randomST, diso, 'disordered')
    randomST_relative = calculate_consurf_difference_psites(randomST_dis, consurf, window_size=5)
    randomST_dis_consurf_references, randomST_dis_consurf = retrieve_ConSurf_score(randomST_dis, consurf)

    disordered_data = [cond_relative, randomST_relative, univ_relative]
    labels = ['Non-phosphorylated S/T', 'Conditional phosphosites', 'Universal phosphosites']

    print(f'conditional: {np.median(cond_relative)}')
    print(f'universal: {np.median(univ_relative)}')
    print(f'Non-phosphorylated S/T: {np.median(randomST_relative)}')

    for db in ['sgd', 'biogrid', 'lanz70', 'lanz90', 'leutert', 'conditional']:
        if db == 'sgd': 
            psites = sgd
        elif db == 'biogrid':
            psites = biogrid
        elif db == 'lanz70':
            psites = lanz70
        elif db == 'lanz90':
            psites = lanz90
        elif db == 'leutert':
            psites = all_psites
        elif db == 'conditional': 
            psites = cond_psites
        print(f'\n{db}')
        all_psites_ST_db = retrieve_references_by_residue_type(psites, sample_residues)
        all_psites_dis_db = retrieve_references_by_order(all_psites_ST_db, diso, 'disordered')
        db_relative = calculate_consurf_difference_psites(all_psites_dis_db, consurf, window_size=5)
        all_dis_consurf_references_db, all_dis_consurf_db = retrieve_ConSurf_score(all_psites_dis_db, consurf)
        print(f'{np.median(db_relative):.3f} of {len(db_relative)}')
        dispensable = estimate_pi_mixture_model(randomST_relative, univ_relative, db_relative)
        print(f'{dispensable:.3f}')

if __name__ == '__main__':
    main()