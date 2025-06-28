'''
ConSurf score distribution among conditional p-sites with number of
perturbations 1, 2-4, 5-7, 8-10. 
Corresponds to Figure 3B
'''

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
    Fig3B = paperDir / 'Figure 3B.jpg'

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
    data_lists = [references_1, references_2_4, references_5_7, references_8_10]
    data_lists = [retrieve_references_by_residue_type(l, sample_residues) for l in data_lists]
    data_lists = [retrieve_references_by_order(l, diso, 'disordered') for l in data_lists]
    data_lists = [retrieve_ConSurf_score(l, consurf)[1] for l in data_lists]
    
    if not Fig3B.is_file():
        plot_consurf_distribution_against_perturbations(data_lists, 
                                                        ['1', '2-4', '5-7', '8-10'], 
                                                        Fig3B, 
                                                        figFmt)
    
if __name__ == '__main__':
    main()