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
    lanz70PKL = paperDir / 'lanz70.pkl'
    PhosphoGIRD_PKL = paperDir / 'PhosphoGIRD.pkl'

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

    print(len(PhosphoGIRD))

    PhosphoGIRD_ST = retrieve_references_by_residue_type(PhosphoGIRD, sample_residues)
    PhosphoGIRD_dis = retrieve_references_by_order(PhosphoGIRD_ST, diso, 'disordered')
    PhosphoGIRD_ord = retrieve_references_by_order(PhosphoGIRD_ST, diso, 'ordered')
    PhosphoGIRD_dis_consurf_references, PhosphoGIRD_dis_consurf = retrieve_ConSurf_score(PhosphoGIRD_dis, consurf)
    PhosphoGIRD_ord_consurf_references, PhosphoGIRD_ord_consurf = retrieve_ConSurf_score(PhosphoGIRD_ord, consurf)

    print(len(PhosphoGIRD_dis_consurf), len(PhosphoGIRD_ord_consurf))

    print(sorted(PhosphoGIRD_dis_consurf))

    print(np.median(PhosphoGIRD_dis_consurf))


if __name__ == '__main__':
    main()