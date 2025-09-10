'''
Overlap between regulated p-sites and universal/conditional p-sites. 
Corresponds to Figure 1C, but does not look elegant, so is finally 
replaced by drawing on PowerPoint. 
'''

import os
import pickle
from pathlib import Path
from plot_tools import plot_reg_overlap
from proteomic_tools import (retrieve_references_by_residue_type, 
							 retrieve_references_by_order)

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
    regPhosPKL = paperDir / 'reguPhos.pkl'

    # output files
    Fig1B = paperDir / 'Figure 1B.jpg'

    phosStres = pickle.load(open(phosStresPKL, 'rb'))
    phosStres = {key: retrieve_references_by_residue_type(references, resType) for key, references in phosStres.items()}
    # phosStres_ord = {key: retrieve_references_by_order(references, diso, 'ordered') for key, references in phosStres.items()}
    # phosStres_dis = {key: retrieve_references_by_order(references, diso, 'disordered') for key, references in phosStres.items()}

    cond = set([reference for key, references in phosStres.items() for reference in references if key < 11])
    univ = set([reference for key, references in phosStres.items() for reference in references if key > 91])
    cond_dis = retrieve_references_by_order(univ, diso, 'disordered')
    univ_dis = retrieve_references_by_order(univ, diso, 'disordered')
    reg = pickle.load(open(regPhosPKL, 'rb'))

    print(len(univ_dis))
    print(len(cond_dis))
    print(len(univ_dis))

if __name__ == '__main__':
    main()
