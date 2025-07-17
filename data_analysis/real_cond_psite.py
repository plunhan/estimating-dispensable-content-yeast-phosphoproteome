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
							 retrieve_ConSurf_score)

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

    print(len(sgd), len(biogrid), len(sgd - biogrid), len(biogrid - sgd))

    sgd = retrieve_references_by_residue_type(sgd, sample_residues)
    sgd = retrieve_references_by_order(sgd, diso, 'disordered')
    biogrid = retrieve_references_by_residue_type(biogrid, sample_residues)
    biogrid = retrieve_references_by_order(biogrid, diso, 'disordered')
    sgd = set(sgd)
    biogrid = set(biogrid)
    print(len(sgd), len(biogrid), len(sgd - biogrid), len(biogrid - sgd))

    '''
    # find a method of binning
    universal = {key: value for key, value in phosStres.items() if key > 91}
    universal = {key: retrieve_references_by_residue_type(value, sample_residues) for key, value in universal.items()}
    universal = {key: retrieve_references_by_order(value, diso, 'disordered') for key, value in universal.items()}
    universal = {key: retrieve_ConSurf_score(value, consurf)[1] for key, value in universal.items()}

    def se_median(data, n_bootstrap=1000, random_state=None):
        rng = np.random.default_rng(random_state)
        medians = np.array([
            np.median(rng.choice(data, size=len(data), replace=True))
            for _ in range(n_bootstrap)
        ])
        return np.std(medians, ddof=1)

    for key, value in universal.items():
        print(f'{key}, median {np.median(value)}, se {se_median(value)}, number of floats {len(value)}')

    l = []
    for key, value in universal.items():
        print(key)
        if key <= 96:
            l += value
        elif 97 <= key <= 100:
            if key == 97:
                print(f'median {np.median(l)}, se {se_median(l)}')
                l = []
            l += value
        elif key == 101:
            print(f'median {np.median(l)}, se {se_median(l)}')
            print(f'median {np.median(value)}, se {se_median(value)}')
    
    '''
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
    data_lists_cond_dis = [retrieve_ConSurf_score(l, consurf)[1] for l in data_lists_cond_dis]
    
    if not Fig4A.is_file():
        plot_consurf_distribution_against_perturbations(data_lists_cond_dis, 
                                                        ['1', '2-4', '5-7', '8-10'], 
                                                        Fig4A, 
                                                        figFmt, 
                                                        (-0.05, 0.5), 
                                                        "Distribution of ConSurf score against number of perturbations\nin disordered regions")

    references_92_96 = []
    references_97_100 = []
    references_101 = []
    for i in range(92, 102): 
        if i == 101:
            references_101.extend(phosStres[i])
        elif 92 <= i <= 96: 
            references_92_96.extend(phosStres[i])
        elif 97 <= i <= 100: 
            references_97_100.extend(phosStres[i])
    data_lists_univ = [references_92_96, references_97_100, references_101]
    data_lists_univ = [retrieve_references_by_residue_type(l, sample_residues) for l in data_lists_univ]
    data_lists_univ_dis = [retrieve_references_by_order(l, diso, 'disordered') for l in data_lists_univ]
    data_lists_univ_dis = [retrieve_ConSurf_score(l, consurf)[1] for l in data_lists_univ_dis]

    if not Fig4B.is_file():
        plot_consurf_distribution_against_perturbations(data_lists_univ_dis, 
                                                        ['92-96', '97-100', '101'], 
                                                        Fig4B, 
                                                        figFmt, 
                                                        (-0.05, 0.5), 
                                                        "Distribution of ConSurf score against number of perturbations\nin disordered regions")
    
if __name__ == '__main__':
    main()