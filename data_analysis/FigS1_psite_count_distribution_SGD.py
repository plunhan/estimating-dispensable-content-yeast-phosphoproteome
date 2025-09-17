'''
ConSurf score comparison between conditional p-sites, universal p-sites, and random S/T. 
'''

import numpy as np
import os
import pandas as pd
import pickle
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path
from plot_tools import plot_consurf_difference
from proteomic_tools import (get_phosphosites_given_perturbations, 
                             retrieve_references_by_residue_type,
                             retrieve_references_by_order,
                             retrieve_ConSurf_score, 
                             sample_random_sites, 
                             calculate_consurf_difference_psites)
from scipy.stats import ranksums

def main(): 

    figFmt = 'jpg'

    sample_residues = 'ST'
    # sample_residues = 'ACDEFGHIKLMNPQRSTVWY'

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
    lanz90PKL = paperDir / 'lanz90.pkl'
    sgdPKL = paperDir / 'SGD.pkl'
    biogridPKL = paperDir / 'BioGRID.pkl'

    # output files
    FigS1 = paperDir / 'Figure S1.jpg'

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
    intermediate_psites = get_phosphosites_given_perturbations(phosStres, list(range(11, 92)))
    all_psites = get_phosphosites_given_perturbations(phosStres, list(range(1, 102)))

    cond_psites = cond_psites.intersection(sgd)
    univ_psites = univ_psites.intersection(sgd)

    cond_psites_ST = retrieve_references_by_residue_type(cond_psites, sample_residues)
    intermediate_psites_ST = retrieve_references_by_residue_type(intermediate_psites, sample_residues)
    univ_psites_ST = retrieve_references_by_residue_type(univ_psites, sample_residues)
    all_psites_ST = retrieve_references_by_residue_type(all_psites, sample_residues)

    cond_psites_dis = retrieve_references_by_order(cond_psites_ST, diso, 'disordered')
    intermediate_psites_dis = retrieve_references_by_order(intermediate_psites_ST, diso, 'disordered')
    univ_psites_dis = retrieve_references_by_order(univ_psites_ST, diso, 'disordered')
    all_psites_dis = retrieve_references_by_order(all_psites_ST, diso, 'disordered')
    cond_psites_ord = retrieve_references_by_order(cond_psites_ST, diso, 'ordered')
    univ_psites_ord = retrieve_references_by_order(univ_psites_ST, diso, 'ordered')
    all_psites_ord = retrieve_references_by_order(all_psites_ST, diso, 'ordered')

    cond_dis_consurf_references, cond_dis_consurf = retrieve_ConSurf_score(cond_psites_dis, consurf)
    intermediate_dis_consurf_references, intermediate_dis_consurf = retrieve_ConSurf_score(intermediate_psites_dis, consurf)
    univ_dis_consurf_references, univ_dis_consurf = retrieve_ConSurf_score(univ_psites_dis, consurf)
    cond_ord_consurf_references, cond_ord_consurf = retrieve_ConSurf_score(cond_psites_ord, consurf)
    univ_ord_consurf_references, univ_ord_consurf = retrieve_ConSurf_score(univ_psites_ord, consurf)
    all_ord_consurf_references, all_ord_consurf = retrieve_ConSurf_score(all_psites_ord, consurf)
    all_dis_consurf_references, all_dis_consurf = retrieve_ConSurf_score(all_psites_dis, consurf)

    exclusions = ultradeep.union(sgd, biogrid, lanz90) # All reported p-sites
    randomST = sample_random_sites(ultradeep, exclusions, sequences, sample_residues)
    randomST_dis = retrieve_references_by_order(randomST, diso, 'disordered')
    randomST_ord = retrieve_references_by_order(randomST, diso, 'ordered')
    randomST_dis_consurf_references, randomST_dis_consurf = retrieve_ConSurf_score(randomST_dis, consurf)
    randomST_ord_consurf_references, randomST_ord_consurf = retrieve_ConSurf_score(randomST_ord, consurf)

    cond_relative = calculate_consurf_difference_psites(cond_psites_dis, consurf, window_size=5)
    univ_relative = calculate_consurf_difference_psites(univ_psites_dis, consurf, window_size=5)
    randomST_relative = calculate_consurf_difference_psites(randomST_dis, consurf, window_size=5)

    ordered_data = [cond_ord_consurf, randomST_ord_consurf, univ_ord_consurf]
    disordered_data = [randomST_relative, cond_relative, univ_relative]
    labels = ['Non-phosphorylated S/T sites', 'Conditional phosphosites', 'Universal phosphosites']

    if not FigS1.is_file():
        plot_consurf_difference(disordered_data, 
                                labels, 
                                FigS1,
                                (-0.5, 0.2), 
                                figFmt, 
                                'disordered',
                                urge_positive=True)

if __name__ == '__main__':
    main()