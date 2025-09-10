'''
ConSurf score comparison between conditional p-sites, universal p-sites, and random S/T. 
'''

import numpy as np
import os
import pandas as pd
import pickle
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path
from plot_tools import plot_consurf_distribution_separate
from proteomic_tools import (get_phosphosites_given_perturbations, 
							 retrieve_references_by_residue_type,
							 retrieve_references_by_order,
							 retrieve_ConSurf_score, 
							 sample_random_sites)

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

    # output files
    Fig3A = paperDir / 'Figure 3A.jpg'
    Fig3B = paperDir / 'Figure 3B.jpg'

    ultradeep = pickle.load(open(ultradeepPKL, 'rb'))
    phosStres = pickle.load(open(phosStresPKL, 'rb'))
    consurf = pickle.load(open(consurfPKL, 'rb'))
    diso = pickle.load(open(disoPKL, 'rb'))
    sequences = pickle.load(open(sequencePKL, 'rb'))
    sgd = pickle.load(open(sgdPKL, 'rb'))
    biogrid = pickle.load(open(biogridPKL, 'rb'))
    lanz90 = pickle.load(open(lanz90PKL, 'rb'))

    # Part I: difference between conditional and universal phosphosites
    cond_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 11)))
    univ_psites = get_phosphosites_given_perturbations(phosStres, list(range(92, 102)))
    all_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 102)))

    for sample_residue in sample_residues:
        cond_psites_ST = retrieve_references_by_residue_type(cond_psites, sample_residue)
        univ_psites_ST = retrieve_references_by_residue_type(univ_psites, sample_residue)
        all_psites_ST = retrieve_references_by_residue_type(all_psites, sample_residues)

        cond_psites_dis = retrieve_references_by_order(cond_psites_ST, diso, 'disordered')
        univ_psites_dis = retrieve_references_by_order(univ_psites_ST, diso, 'disordered')
        all_psites_dis = retrieve_references_by_order(all_psites_ST, diso, 'disordered')

        cond_dis_consurf_references, cond_dis_consurf = retrieve_ConSurf_score(cond_psites_dis, consurf)
        univ_dis_consurf_references, univ_dis_consurf = retrieve_ConSurf_score(univ_psites_dis, consurf)

        print(sample_residue)
        print(f'Conditional: {len(cond_dis_consurf)}')
        print(f'Universal: {len(univ_dis_consurf)}')

        exclusions = ultradeep.union(sgd, biogrid, lanz90) # All reported p-sites
        randomST = sample_random_sites(ultradeep, exclusions, sequences, sample_residue)
        randomST_dis = retrieve_references_by_order(randomST, diso, 'disordered')
        randomST_dis_consurf_references, randomST_dis_consurf = retrieve_ConSurf_score(randomST_dis, consurf)

        disordered_data = [cond_dis_consurf, randomST_dis_consurf, univ_dis_consurf]
        labels = ['Conditional phosphosites', 'Non-phosphorylated S/T', 'Universal phosphosites']

        if sample_residue == 'S' and not Fig3A.is_file():
            plot_consurf_distribution_separate(disordered_data, 
                                               labels, 
                                               Fig3A, 
                                               figFmt, 
                                               'Distribution of ConSurf score for serines',
                                               sample_residue,)

        if sample_residue == 'T' and not Fig3B.is_file():
            plot_consurf_distribution_separate(disordered_data, 
                                               labels, 
                                               Fig3B, 
                                               figFmt, 
                                               'Distribution of ConSurf score for threonines',
                                               sample_residue)

if __name__ == '__main__':
    main()