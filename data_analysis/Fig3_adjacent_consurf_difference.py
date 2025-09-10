'''
Calculate the difference of ConSurf score between a residue and its adjacent residues (+2, +1, -1, and -2)

Comparison of evolutionary rate between:
   Conditional phosphosites
   Universal phosphosites
   Random S/T
   Random S/T that never been reported by any phosphoproteome
'''

import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from plot_tools import plot_consurf_difference
from proteomic_tools import (get_phosphosites_given_perturbations, 
                             retrieve_references_by_order,
                             retrieve_references_by_residue_type,
                             sample_random_sites, 
                             calculate_consurf_difference_psites)
from scipy.stats import ranksums, ttest_ind

def main():

    sample_residues = set('ST')

    figFmt = 'jpg'

    window_size = 5

    dataDir = Path('../../data')

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
    randomST_dis_consurf_distribution = paperDir / 'randomST_dis_consurf_distribution.jpg'
    Fig3 = paperDir / 'Figure 3.jpg'
    FigS2B = paperDir / 'Figure S2B.jpg'

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

    cond_psites_ST = retrieve_references_by_residue_type(cond_psites, sample_residues)
    univ_psites_ST = retrieve_references_by_residue_type(univ_psites, sample_residues)
    all_psites_ST = retrieve_references_by_residue_type(all_psites, sample_residues)

    cond_psites_dis = retrieve_references_by_order(cond_psites_ST, diso, 'disordered')
    univ_psites_dis = retrieve_references_by_order(univ_psites_ST, diso, 'disordered')
    all_psites_dis = retrieve_references_by_order(all_psites_ST, diso, 'disordered')

    cond_diff_dis = calculate_consurf_difference_psites(cond_psites_dis, consurf, window_size)
    univ_diff_dis = calculate_consurf_difference_psites(univ_psites_dis, consurf, window_size)

    print(f'Conditional: {np.median(cond_diff_dis)}')
    print(f'Universal: {np.median(univ_diff_dis)}')
    _, pvalue = ranksums(cond_diff_dis, univ_diff_dis)
    print(pvalue)

    
    # Ordered regions are not considered, but here is the code. 
    cond_psites_ord = retrieve_references_by_order(cond_psites_ST, diso, 'ordered')
    univ_psites_ord = retrieve_references_by_order(univ_psites_ST, diso, 'ordered')
    all_psites_ord = retrieve_references_by_order(all_psites_ST, diso, 'ordered')

    cond_diff_ord = calculate_consurf_difference_psites(cond_psites_ord, consurf, window_size)
    univ_diff_ord = calculate_consurf_difference_psites(univ_psites_ord, consurf, window_size)

    print(f'Conditional: {np.median(cond_diff_ord)}')
    print(f'Universal: {np.median(univ_diff_ord)}')
    _, pvalue = ranksums(cond_diff_ord, univ_diff_ord)
    print(pvalue)
    

    exclusions = ultradeep.union(sgd, biogrid, lanz90) # All reported p-sites
    randomST = sample_random_sites(ultradeep, exclusions, sequences, sample_residues)
    randomST_dis = retrieve_references_by_order(randomST, diso, 'disordered')
    randomST_ord = retrieve_references_by_order(randomST, diso, 'ordered')
    randomST_diff_dis = calculate_consurf_difference_psites(randomST_dis, consurf, window_size)
    randomST_diff_ord = calculate_consurf_difference_psites(randomST_ord, consurf, window_size)
    print(np.median(randomST_diff_ord))

    ordered_data = [randomST_diff_ord, cond_diff_ord, univ_diff_ord]
    disordered_data = [randomST_diff_dis, cond_diff_dis, univ_diff_dis]
    labels = ['Non-phosphorylated S/T', 'Conditional phosphosites', 'Universal phosphosites']

    print(ranksums(randomST_diff_ord, cond_diff_ord))
    print(ranksums(univ_diff_ord, randomST_diff_ord))
    print(ranksums(univ_diff_ord, cond_diff_ord))

    if not Fig3.is_file():
        plot_consurf_difference(disordered_data, 
                                labels, 
                                Fig3, 
                                (-0.4, 0.15),
                                figFmt)

    if not FigS2B.is_file():
        plot_consurf_difference(ordered_data, 
                                labels, 
                                FigS2B, 
                                (-0.15, 0.15),
                                figFmt)
    
if __name__ == '__main__':
    main()