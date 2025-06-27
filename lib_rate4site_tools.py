import os
import pandas as pd
import subprocess
from pathlib import Path
from tqdm import tqdm
from typing import Union
from Bio import SeqIO

def prepare_rate4site_files(proteome_path: Union[str, Path], 
                            rate4siteDir: Union[str, Path], 
                            alignmentDir: Union[str, Path]) -> list[str]:
    '''
    Prepare input files for rate4site. 
    
    Args:
        proteome_path (Union[str, Path]): Path to S. cerevisiae ORF fasta file. 
        rate4siteDir (Union[str, Path]): Directory to contain input files for rate4site. 
        alignmentDir (Union[str, Path]): Multiple-sequence alignment directory. 

    Returns:
        list[str]: List of S. cerevisiae ORFs that have an ortholog in all closely-related yeast species. 
    '''
    orfs = []
    with open(proteome_path, 'r') as file: 
        for line in file: 
            if line.startswith('>'): 
                protein = line.split('>')[1].split()[0]
                orfs.append(protein)
    output_orfs = []
    for orf in tqdm(orfs): 
        alignmentFile = alignmentDir / orf / (orf + '.aa.aln')
        if alignmentFile.is_file(): 
            output_orfs.append(orf)
            orfDir = rate4siteDir / orf
            if not orfDir.exists():
                os.mkdir(orfDir)
            alignmentFile_copy = orfDir / (orf + '.aa.aln')
            if not alignmentFile_copy.is_file(): 
                with open(alignmentFile, 'r') as original, open(alignmentFile_copy, 'w') as corrected:
                    records = SeqIO.parse(original, 'fasta')
                    for record in records:
                        record.id = record.description.split()[1][:4].capitalize()
                        record.description = record.id
                        SeqIO.write(record, corrected, 'fasta')
    return output_orfs

def run_rate4site(input_ls: list[str], 
                  tree_file: Union[str, Path], 
                  inDir: Union[str, Path], 
                  outDir: Union[str, Path], 
                  rate4sitePath: Union[str, Path]) -> None: 
    '''
    Run rate4site on each S. cerevisiae ORF. 

    Args:
        input_ls (list[str]): S. cerevisiae ORFs of interest. 
        tree_file (Union[str, Path]): A specified tree. 
        inDir (Union[str, Path]): Parent directory of alignment results. 
        outDIr (Union[str, Path]): Parent directory of ConSurf results. 
        rate4sitePath (Union[str, Path]): Path to executable rate4site software. 

    Returns:
        None
    '''
    print("-------- Run rate4Site -------- ")
    if not outDir.exists():
        os.makedirs(outDir)
    for orf in tqdm(input_ls): 
        aln_file = inDir / orf / (orf + '.aa.aln')
        outfile = outDir / (orf + '.consurf.txt')
        if aln_file.is_file() and (not outfile.is_file()):
            cmd = str(rate4sitePath) + ' -s ' + str(aln_file) + ' -t ' + str(tree_file) + ' -o ' + str(outfile)
            # cmd = [str(rate4sitePath), '-s', str(aln_file), '-t', str(tree_file), '-o', str(outfile)]
            subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,cwd=outDir)