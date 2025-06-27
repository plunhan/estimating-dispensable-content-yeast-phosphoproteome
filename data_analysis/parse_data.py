# Parse the following results into .pkl files:
#   Ultradeep phosphoproteome
#   Perturbation-specific phosphoproteome
#   SGD phosphoproteome
#   BioGRID phosphoproteome
#   S. cerevisiae sequences
#   Rate4Site results (ConSurf scores)
#   DISOPRED results (ordered/disorderd regions)

import os
import pickle
from pathlib import Path

from lib_import_tools import map_protein_id_to_locus_id, map_gene_name_to_locus_id
from proteomic_tools import (parse_ultradeep, 
                             parse_stress_phosphosites, 
                             parse_disopred,
                             parse_consurf,
                             parse_fasta,
                             parse_sgd,
                             parse_biogrid,
                             parse_lanz,
                             parse_functional_phosphosites, 
                             parse_regulated_phosphosites)

def main():

    dataDir = Path('../data')

    extDir = dataDir / 'external'

    procDir = dataDir / 'processed'

    figDir = procDir / 'figures'

    LeutertDir = extDir / 'Leutert'

    disoDir = procDir / 'all_disopred'

    consurfDir = procDir / 'ConSurf'

    paperDir = procDir / 'paper'

    if not paperDir.exists():
        os.mkdir(paperDir)

    # input files
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    ultradeepFile = LeutertDir / 'ultradeep_reference_phosphoproteome.xlsx'
    phosStresFile = LeutertDir / 'quantitative_phosphosites.csv'
    sequenceFile = extDir / 'Scer.aa'
    sgdFile = extDir / 'alliancemine_results_2025-03-26T11-40-38.tsv'
    biogridFile = extDir / 'BIOGRID-PTM-4.4.235.ptmtab.txt'
    funcPhosFile = extDir / '41587_2021_1051_MOESM9_ESM.xlsx'
    lanzFile = extDir / 'embr202051121-phosphoproteome.xlsx'
    diffFile = LeutertDir / 'Differential expression all.xlsx'

    # output files
    ultradeepPKL = paperDir / 'ultradeep_reference_phosphoproteome.pkl'
    phosStresPKL = paperDir / 'quantitative_phosphosites.pkl'
    consurfPKL = paperDir / 'consurf_all.pkl'
    disoPKL = paperDir / 'diso_all.pkl'
    sequencePKL = paperDir / 'Scer_seq.pkl'
    sgdPKL = paperDir / 'SGD.pkl'
    biogridPKL = paperDir / 'BioGRID.pkl'
    funcPhosPKL = paperDir / 'funcPhos.pkl'
    lanzPKL = paperDir / 'lanz.pkl'
    regPhosPKL = paperDir / 'reguPhos.pkl'

    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    IDMappingDict2 = map_gene_name_to_locus_id(IDMappingFile)

    if not ultradeepPKL.is_file():
        ultradeep = parse_ultradeep(ultradeepFile)
        pickle.dump(ultradeep, open(ultradeepPKL, 'wb'))

    if not phosStresPKL.is_file():
        phosStres = parse_stress_phosphosites(phosStresFile)
        pickle.dump(phosStres, open(phosStresPKL, 'wb'))

    if not disoPKL.is_file():
        diso = parse_disopred(disoDir, IDMappingDict)
        pickle.dump(diso, open(disoPKL, 'wb'))

    if not consurfPKL.is_file():
        consurf = parse_consurf(consurfDir)
        pickle.dump(consurf, open(consurfPKL, 'wb'))

    if not sequencePKL.is_file():
        sequences = parse_fasta(sequenceFile)
        pickle.dump(sequences, open(sequencePKL, 'wb'))

    if not sgdPKL.is_file():
        sgd = parse_sgd(sgdFile)
        pickle.dump(sgd, open(sgdPKL, 'wb'))

    if not biogridPKL.is_file():
        biogrid = parse_biogrid(biogridFile)
        pickle.dump(biogrid, open(biogridPKL, 'wb'))

    if not funcPhosPKL.is_file():
        funcPhos = parse_functional_phosphosites(funcPhosFile, IDMappingDict2)
        pickle.dump(funcPhos, open(funcPhosPKL, 'wb'))

    if not lanzPKL.is_file():
        if not sequencePKL.is_file():
            sequences = parse_fasta(sequenceFile)
            pickle.dump(sequences, open(sequencePKL, 'wb'))
        sequences = pickle.load(open(sequencePKL, 'rb'))
        lanz = parse_lanz(lanzFile, sequences)
        pickle.dump(lanz, open(lanzPKL, 'wb'))

    if not regPhosPKL.is_file():
        regPhos = parse_regulated_phosphosites(diffFile)
        pickle.dump(regPhos, open(regPhosPKL, 'wb'))

if __name__ == '__main__':
    main()