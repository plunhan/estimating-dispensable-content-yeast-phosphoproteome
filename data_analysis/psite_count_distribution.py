'''
Plot distribution of count of p-sites against the number of perturbations in which 
they are detected. 
Correspond to Figures 1A and 1B
'''

import os
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ranksums, ttest_ind

from plot_tools import plot_count_per_perturbation
from proteomic_tools import retrieve_references_by_residue_type, retrieve_references_by_order

def main():

    resType = 'ST'

    figFmt = 'jpg'

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

    # output files
    Fig1A = paperDir / 'Figure 1A.jpg'
    Fig1B = paperDir / 'Figure 1B.jpg'

    diso = pickle.load(open(disoPKL, 'rb'))
    phosStres = pickle.load(open(phosStresPKL, 'rb'))
    phosStres = {key: retrieve_references_by_residue_type(references, resType) for key, references in phosStres.items()}
    phosStres_ord = {key: retrieve_references_by_order(references, diso, 'ordered') for key, references in phosStres.items()}
    phosStres_dis = {key: retrieve_references_by_order(references, diso, 'disordered') for key, references in phosStres.items()}

    if not Fig1A.is_file(): 
        plot_count_per_perturbation(phosStres_ord, 
                                    Fig1A, 
                                    'Detection frequency of phosphosites in ordered\nregions across perturbations',
                                    figFmt)
    if not Fig1B.is_file(): 
        plot_count_per_perturbation(phosStres_dis, 
                                    Fig1B, 
                                    'Detection frequency of phosphosites in disordered\nregions across perturbations',
                                    figFmt)

if __name__ == '__main__':
    main()
