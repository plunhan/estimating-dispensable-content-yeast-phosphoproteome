'''
Pairwise sequence alignment between each protein of S. cerevisiae 
and each protein of a closely-related yeast species. 
'''
import os
from Bio import SeqIO
from lib_blast_tools import process_uniprot_fasta, run_blast, coverage_cutoff_on_BLAST
from lib_import_tools import map_protein_id_to_locus_id
from pathlib import Path

def main():

    coverage_cutoff = 0.5
    evalue_cutoff = 1e-5
    threads = 8

    dataDir = Path('../data')

    extDir = dataDir / 'external'

    procDir = dataDir / 'processed'

    if not procDir.exists(): 
        os.mkdir(procDir)

    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)

    # process uniprot fasta file
    uProt = extDir / 'UP000002311_559292.fasta'
    proProt = extDir / 'Scer.aa'
    uDNA = extDir / 'UP000002311_559292_DNA.fasta'
    proDNA = extDir / 'Scer.nt'

    process_uniprot_fasta(uProt, IDMappingDict, proProt)
    process_uniprot_fasta(uDNA, IDMappingDict, proDNA)
    
    species_ls = ['Sbay', 'Smik', 'Spar', 'Ncas', 'Cgla', 'Agos', 'Klac', 'Calb']
    for species in species_ls: 
        databaseName = 'Scer.aa'
        outPath1 = procDir / ('Scer-' + species + '-blast_stats.best.txt')
        outPath2 = procDir / ('Scer-' + species + '-blast_stats_coverageCutoff.best.txt')
        queryFile = extDir / (species + '.aa')
        if not outPath1.is_file(): 
            print('running blast')
            run_blast(queryFile,
                      databaseName,
                      extDir,
                      procDir,
                      outPath1,
                      threads = threads,
                      evalue_cutoff = evalue_cutoff)
        if not outPath2.is_file(): 
            coverage_cutoff_on_BLAST(outPath1, outPath2, coverage_cutoff=coverage_cutoff)

if __name__ == '__main__':
    main()