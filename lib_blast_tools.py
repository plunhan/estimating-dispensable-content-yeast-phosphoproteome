import os
import subprocess
from pathlib import Path
from typing import Union

def process_uniprot_fasta(inPath: Union[str, Path], 
                          IDMappingDict: dict[str, str], 
                          outPath: Union[str, Path]): 
    '''
    Replace UniProt IDs with Systematic IDs in a fasta file. 

    Args:
        inPath (Union[str, Path]): Path to fasta with UniProt IDs in header. 
        IDMappingDict (dict[str, str]): Dictionary mapping UniProt IDs to Systematic IDs. 
        outPath (Union[str, Path]): Path to fasta with Systematic IDs in header. 

    Returns:
        None
    '''
    filtered_sequences = []
    for record in SeqIO.parse(inPath, 'fasta'):
        original_header = record.id
        id_part = original_header.split('|')[1]
        if id_part in IDMappingDict.keys():
            record.id = IDMappingDict[id_part]
            filtered_sequences.append(record)
    SeqIO.write(filtered_sequences, outPath, 'fasta')

def run_blast(queryFile: Union[str, Path],
              databaseName: str, 
              extDir: Union[str, Path], 
              procDir: Union[str, Path],
              outPath: Union[str, Path],
              threads = 8, 
              evalue_cutoff = 1e-5) -> None:
    '''
    Run BLAST sequence alignment using a query file against a database. 

    Args: 
        queryFile (Union[str, Path]): Query sequence file. 
        databaseName (str): The name of database that BLAST aligns against. 
        extDir (Union[str, Path]): Directory of fasta file for building the database. 
        procDir (Union[str, Path]): Directory of database. 
        outPath (Union[str, Path]): Directory of the results. 
        threads (int): Number of GPUs for computing. 
        evalue_cutoff (float): E-value cutoff, worse than which the alignment is discarded. 

    Returns:
        None
    '''
    databaseFile = extDir / databaseName
    databaseResult = procDir / (databaseName + '.pin')
    if not databaseResult.is_file(): 
        cmd_makedb = ['makeblastdb', '-in', databaseFile, '-dbtype', 'prot', '-out', str(procDir / databaseName)]
        subprocess.run(cmd_makedb)
    if not outPath.is_file(): 
        fmt = '"6 sseqid qseqid slen qlen length sstart send qstart qend nident positive gaps evalue"'
        os.system('blastp -query ' + str(queryFile) +
              ' -db ' + str(procDir / databaseName) +
              ' -evalue ' + str(evalue_cutoff) +
              ' -outfmt ' + fmt +
              ' -num_threads ' + str(threads) + ' > ' + str(outPath) +
              ' 2>/dev/null')
        # # The following commands did not word for my computer so I used os.system instead
        # cmd_align = ['blastp', '-query', str(queryFile), '-db', 
        #              str(procDir / databaseName), '-evalue', str(evalue_cutoff), 
        #              '-outfmt', fmt, '-num_threads', str(threads), 
        #              '>', str(outPath), '2>/dev/null']
        # subprocess.run(cmd_align)

def coverage_cutoff_on_BLAST(inPath: Union[str, Path], 
                             outPath: Union[str, Path], 
                             coverage_cutoff = 0.5) -> None: 
    '''
    Filter BLAST results based on alignment coverage.

    Args:
        inPath (Union[str, Path]): Path to the input file containing BLAST results.
        outPath (Union[str, Path]): Path to the output file for filtered results.
        coverage_cutoff (float): Minimum alignment coverage required for both query and subject (default is 0.5, i.e., 50%).

    Returns:
        None
    '''
    with open(inPath, 'r') as fin, open(outPath, 'w') as fout:
        for line in fin.readlines():
            line = line.split('\t')
            target = line[0]
            query = line[1]
            qlen = float(line[2])
            slen = float(line[3])
            qstart = float(line[5])
            qend = float(line[6])
            sstart = float(line[7])
            send = float(line[8])
            positives = int(line[10])
            gaps = int(line[11])
            evalue = float(line[12])
            if ( ( (qend - qstart + 1 - gaps) / qlen > coverage_cutoff ) and \
                 ( (send - sstart + 1 - gaps) / slen > coverage_cutoff )): 
                fout.write(target + '\t' + query + '\t' + str(evalue) + '\n')