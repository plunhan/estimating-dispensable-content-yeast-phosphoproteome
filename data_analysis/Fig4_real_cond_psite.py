'''
ConSurf score distribution among conditional p-sites with number of
perturbations 1, 2-4, 5-7, 8-10. 
Corresponds to Figure 3B
'''

import numpy as np
import os
import pandas as pd
import pickle
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path
from plot_tools import plot_consurf_distribution, plot_consurf_distribution_against_perturbations
from proteomic_tools import (retrieve_references_by_residue_type,
							 retrieve_references_by_order,
							 retrieve_ConSurf_score, 
                             calculate_consurf_difference_psites)

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

    # output files
    Fig4A = paperDir / 'Figure 4A.jpg'
    Fig4B = paperDir / 'Figure 4B.jpg'

    ultradeep = pickle.load(open(ultradeepPKL, 'rb'))
    phosStres = pickle.load(open(phosStresPKL, 'rb'))
    consurf = pickle.load(open(consurfPKL, 'rb'))
    diso = pickle.load(open(disoPKL, 'rb'))
    sequences = pickle.load(open(sequencePKL, 'rb'))
    sgd = pickle.load(open(sgdPKL, 'rb'))
    biogrid = pickle.load(open(biogridPKL, 'rb'))

    references_1 = []
    references_2_4 = []
    references_5_7 = []
    references_8_10 = []
    for i in range(1, 11): 
        if i == 1:
            references_1.extend(phosStres[i])
        elif 1 < i <= 4: 
            references_2_4.extend(phosStres[i])
        elif 4 < i <= 7: 
            references_5_7.extend(phosStres[i])
        elif i > 7: 
            references_8_10.extend(phosStres[i])
    data_lists_cond = [references_1, references_2_4, references_5_7, references_8_10]
    data_lists_cond = [retrieve_references_by_residue_type(l, sample_residues) for l in data_lists_cond]
    data_lists_cond_dis = [retrieve_references_by_order(l, diso, 'disordered') for l in data_lists_cond]
    data_lists_cond_dis = [calculate_consurf_difference_psites(l, consurf, window_size=5) for l in data_lists_cond_dis]

    
    if not Fig4A.is_file():
        plot_consurf_distribution_against_perturbations(data_lists_cond_dis, 
                                                        ['1', '2-4', '5-7', '8-10'], 
                                                        Fig4A, 
                                                        figFmt, 
                                                        (-0.5, 0.05), 
                                                        "Distribution of relative ConSurf score against number of perturbations\nin disordered regions")

    references_92_94 = []
    references_95_97 = []
    references_98_100 = []
    references_101 = []
    for i in range(92, 102): 
        if i == 101:
            references_101.extend(phosStres[i])
        elif 92 <= i <= 94: 
            references_92_94.extend(phosStres[i])
        elif 95 <= i <= 97: 
            references_95_97.extend(phosStres[i])
        elif 98 <= i <= 100: 
            references_98_100.extend(phosStres[i])
    data_lists_univ = [references_92_94, references_95_97, references_98_100, references_101]
    data_lists_univ = [retrieve_references_by_residue_type(l, sample_residues) for l in data_lists_univ]
    data_lists_univ_dis = [retrieve_references_by_order(l, diso, 'disordered') for l in data_lists_univ]
    data_lists_univ_dis = [calculate_consurf_difference_psites(l, consurf, window_size=5) for l in data_lists_univ_dis]

    print([np.median(l) for l in data_lists_univ_dis])

    if not Fig4B.is_file():
        plot_consurf_distribution_against_perturbations(data_lists_univ_dis, 
                                                        ['92-94', '95-97', '98-100', '101'], 
                                                        Fig4B, 
                                                        figFmt, 
                                                        (-0.5, 0.05), 
                                                        "Distribution of relative ConSurf score against number of perturbations\nin disordered regions")
    
if __name__ == '__main__':
    main()