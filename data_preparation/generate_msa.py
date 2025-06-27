'''
Generate multiple-sequence alignment using ClustalO on BLAST results. 
'''
import os
from pathlib import Path
from lib_msa_tools import generate_yeast_orthologs_all, build_yeast_alignment

def main(): 

    dataDir = Path('../data')

    extDir = dataDir / 'external'

    procDir = dataDir / 'processed'

    clustalwPath = Path('../clustalo/clustal-omega-1.2.3-macosx')

    outPath = procDir / 'codon_alignment.txt'

    orthologs = generate_yeast_orthologs_all(procDir)

    if not outPath.is_file(): 
        build_yeast_alignment(extDir, procDir, orthologs, outPath, clustalwPath)

if __name__ == '__main__':
    main()