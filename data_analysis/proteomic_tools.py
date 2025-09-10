import os
import random
from collections import Counter
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

def parse_ultradeep(ultradeep: Union[str, Path]) -> set[str]:
    '''
    Parse Ultradeep phosphoproteome and return the set of all phosphosites. 

    Args: 
        ultradeep (Union[str, Path]): Path to the Ultradeep phosphoproteome file. 

    Returns:
        set[str]: All phosphosite references in the Ultradeep phosphoproteome. 
    '''
    phos = pd.read_excel(ultradeep, sheet_name='deep_phospho_psites')
    return set(phos['reference'])

def parse_stress_phosphosites(inPath: Union[str, Path]) -> dict[int, set[str]]: 
    '''
    Parse perturbation-specific phosphoproteome. 

    Args:
        inPath (Union[str, Path]): Path to the perturbation-specific phosphoproteome file.

    Returns:
        dict[int, set[str]]: A dictionary mapping the number of perturbations to phosphosite references.
    '''
    df = pd.read_csv(inPath)
    result_d = {i: set() for i in range(1, 102)}
    record_d = {}
    for index, row in df.iterrows():
        if row['reference'] not in record_d:
            record_d[row['reference']] = set()
        record_d[row['reference']].add(row['treatment_id'])
    for reference, treatment in record_d.items():
        result_d[len(treatment)].add(reference)
    return result_d

def get_phosphosites_given_perturbations(phosStres: dict[int, set[str]], 
                                         pert_ls: list[int]) -> set[str]:
    '''
    Return phosphosite references given the numbers of perturbations. 

    Args:
        phosStres (dict[int, set[str]]): A dictionary mapping the number of perturbations and the set of phosphosite references
        pert_ls (list[int]): A list of the numbers of perturbations. 

    Returns:
        set[str]: Phosphosite references given the numbers of perturbations. 
    '''
    references = set()
    for pert in pert_ls:
        if pert not in phosStres: 
            print(f'{pert} is not in the dictionary.')
        references = references.union(phosStres[pert])
    return references

def parse_sgd(inPath: Union[str, Path]) -> set[str]:
    '''
    Return phosphosite references in SGD. 

    Args:
        inPath (Union[str, Path]): Path to the SGD phosphoproteome. 

    Returns:
        set[str]: SGD Phosphosite references. 
    '''
    sgd = pd.read_table(inPath)
    sgd.columns = ['SGDID', 'Systematic ID', 'Length', 'Site', 'Type', 'Modifier', 'PubMedID']
    sgd = sgd[(sgd['Type'] == 'phosphorylated residue')]
    sgd['reference'] = sgd[['Systematic ID', 'Site']].apply('_'.join, axis=1)
    return set(sgd['reference'])

def parse_biogrid(inPath: Union[str, Path]) -> set[str]:
    '''
    Return phosphosite references in BioGRID. 

    Args:
        inPath (Union[str, Path]): Path to the BioGRID phosphoproteome. 

    Returns:
        set[str]: BioGRID Phosphosite references. 
    '''
    biogrid = pd.read_table(inPath)
    biogrid = biogrid[(biogrid['Organism ID'] == 559292) &
                      (biogrid['Post Translational Modification'] == 'Phosphorylation')]
    biogrid['Position'] = biogrid['Position'].astype(str)
    biogrid['site'] = biogrid[['Residue', 'Position']].apply(''.join, axis=1)
    biogrid['reference'] = biogrid[['Systematic Name', 'site']].apply('_'.join, axis=1)
    return set(biogrid['reference'])

def parse_disopred(disorder_dir: Union[str, Path], 
                   IDMappingDict: dict[str, str]) -> dict[str, dict[int, str]]:
    '''
    Parse DISOPRED outputs. 

    Args:
        disorder_dir (Union[str, Path]): Path to the folder containing DISOPRED outputs. 
        IDMappingDict (dict): Dictionary mapping UniProt IDs to Locus IDs. 

    Returns:
        dict[str, dict[int, str]]: A dictionary mapping ORFs to a dictionary mapping position to disordered or ordered regions.
    '''
    return_dict = dict()
    for fname in os.listdir(disorder_dir):
        if not fname.endswith(".diso"):
            continue
        protein = fname.split(".diso")[0]
        if protein not in IDMappingDict: 
            continue
        orf = IDMappingDict[protein]
        return_dict[orf] = {}
        with open(os.path.join(disorder_dir, fname), 'r') as fin: 
            for line in fin: 
                if line.startswith('#'): 
                    continue
                columns = line.split()
                position = int(columns[0])
                amino_acid = columns[1]
                marker = columns[2]
                if marker == "*": 
                    return_dict[orf][position] = "disordered" 
                elif marker == '.': 
                    return_dict[orf][position] = "ordered"
    return return_dict

def parse_fasta(fasta_file: Union[str, Path]) -> dict[str, str]:
    '''
    Parse a FASTA file and return a dictionary mapping protein ID to sequence.

    Args:
        fasta_file (Union[str, Path]): Path to the input FASTA file. 

    Returns:
        dict[str, str]: A dictionary mapping protein ID to sequence. 
    '''
    sequences = {}
    with open(fasta_file) as f:
        curr_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_id is not None:
                    sequences[curr_id] = "".join(seq_lines)
                # Take the first word after ">"
                curr_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if curr_id is not None:
            sequences[curr_id] = "".join(seq_lines)
    return sequences

def parse_consurf(consurf_folder: Union[str, Path]) -> dict[str, dict[int, float]]:
    '''
    Parse ConSurf score files. Lines starting with "#" are skipped.
    
    Args:
        consurf_folder (Union[str, Path]): Path to the folder containing all ConSurf results calculated by Rate4Site. 

    Returns:
        dict[str, dict[int, float]]: A dictionary mapping protein ID to a dictionary mapping position to the ConSurf score.
    '''
    consurf = {}
    for fname in os.listdir(consurf_folder):
        if not fname.endswith(".consurf.txt"):
            continue
        orf = fname.split(".consurf.txt")[0]
        consurf[orf] = {}
        with open(os.path.join(consurf_folder, fname), 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                pos = int(parts[0])
                score = float(parts[2])
                consurf[orf][pos] = score
    return consurf

def parse_regulated_phosphosites(inPath: Union[str, Path]) -> set[str]:
    '''
    Return phosphosite references with a >2-fold change signaling responses, either upregulated or downregulated, with adjusted p-value < 0.05.

    Args:
        inPath (Union[str, Path]): Path to the input file. 

    Returns:
        set[str]: Phosphosite references of regulated phosphosites. 
    '''
    dif = pd.read_excel(inPath, sheet_name='p_site_diff_reg')
    regulated = dif[(dif.adj_p_value < 0.05) & (dif.fc_log2.abs() > 2)]
    return set(regulated['reference'])

def parse_functional_phosphosites(funcPhos: Union[str, Path], 
                                  IDMappingDict: dict[str, str]) -> set[str]:
    '''
    Return phosphosite references that are thought to be "functional". 

    Args:
        funcPhos (Union[str, Path]): Path to functional phosphosite file. 
        IDMappingDict (dict[str, str]): Dictionary mapping yeast protein name to locus ID. 

    Returns:
        set[str]: Phosphosite references after protein ids are mapped. 
    '''
    df = pd.read_excel(funcPhos, sheet_name='go_enrichment_fil_table1')
    references = df['Phosphomutant'].values
    result = set()
    for reference in references:
        protein_site = reference.split('_')[0]
        print(protein_site)
        protein, site = protein_site.rsplit('-', 1)
        if protein not in IDMappingDict:
            if protein in IDMappingDict.values():
                locus = protein
            else:
                continue
        else:
            locus = IDMappingDict[protein]
        site = site[:-1]
        reference_mapped = '_'.join([locus, site])
        result.add(reference_mapped)
    return result

def parse_lanz(lanzFile: Union[str, Path], 
               sequences: dict[str, str], 
               quality: set[str]) -> set[str]:
    '''
    Parse lanz's phosphoproteome (Source: Lanz et al.)

    Args:
        lanzFile (Union[str, Path]): path to lanz's phosphoproteome dataset. 
        sequences (dict[str, str]): dictionary mapping protein name to sequence. 

    Returns:
        set[str]: phosphosite references for in lanz's phosphoproteome dataset. 
    '''
    df = pd.read_excel(lanzFile, sheet_name='Aggregate phosphoproteome')
    df = df[df['Phosphosite Localization'].isin(quality)]
    references_unmapped = df['Gene_Site (unique identifier)']
    references_mapped = set()
    for reference in references_unmapped:
        protein, pos = reference.split('_')
        pos = int(pos)
        if protein not in sequences or pos > len(sequences[protein]) or sequences[protein][pos-1] not in 'STY':
            continue
        aa = sequences[protein][pos-1]
        references_mapped.add(f'{protein}_{aa}{pos}')
    print(f"There are {len(references_mapped)} in Lanz's phosphoproteome.")
    return references_mapped

def parse_phosphogrid_functional(PhosphoGIRDFile: Union[str, Path], 
                                 IDMappingDict: dict[str, str]) -> set[str]:
    df = pd.read_excel(PhosphoGIRDFile, skiprows=8)
    df = df[(df['Function'] == 'Yes') & (df['S,T,Y'] != 'Y')]
    df['Locus'] = df['SGD Gene'].map(IDMappingDict)
    df['reference'] = df.apply(lambda row: f"{row['Locus']}_{row['S,T,Y']}{row['Residue #']}", axis=1)
    return set(df['reference'])

def retrieve_references_by_residue_type(references: set[str],
                                        aatype: set[str]) -> set[str]:
    '''
    Return phosphosite references on specific amino acid type(s).

    Args:
        references (set[str]): Phosphosite references to be filtered. 
        aatype (set[str]): A set of amino acid type to keep. 

    Returns:
        set[str]: Phosphosite references after filtration.
    '''
    references_filtered = set()
    for reference in references: 
        protein, site = reference.split('_')
        aa = site[0]
        if aa in aatype: 
            references_filtered.add(reference)
    return references_filtered

def retrieve_references_by_given_proteins(references: set[str], 
                                          protein_set: set[str]) -> set[str]:
    '''
    Return phosphosite references in specific set of proteins.

    Args:
        references (set[str]): Set of phosphosite references to be filtered. 
        protein_set (set[str]): Set of proteins in which phosphosites are retained. 

    Returns:
        set[str]: A set of phosphosite references after filtration. 
    '''
    result = set()
    for reference in references:
        protein, site = reference.split('_')
        if protein in protein_set: 
            result.add(reference)
    return result

def retrieve_references_by_order(references: set[str], 
                                 order_d: dict[str, dict[int, str]], 
                                 order_type: str) -> set[str]:
    '''
    Return phosphosite references within specific regions (ordered/disordered).

    Args:
        references (set[str]): Set of phosphosite references to be filtered. 
        order_d (dict[str, dict[int, str]]): A dictionary mapping ORFs to a dictionary mapping position to disordered or ordered regions.
        order_type (str): A string indicating which region to retain. Can only be either "ordered" or "disordered". 

    Returns:
        set[str]: A set of phosphosite references after filtration. 
    '''
    result = set()
    for reference in references:
        protein, site = reference.split('_')
        aa = site[0]
        pos = int(site[1:])
        if protein not in order_d or pos not in order_d[protein]:
            continue
        if order_d[protein][pos] == order_type:
            result.add(reference)
    return result

def retrieve_ConSurf_score(references: set[str], 
                           consurf_d: dict[str, dict[int, float]]) -> list[list[float]]:
    '''
    Return ConSurf scores for a given set of phosphosite references.

    Args:
        references (set[str]): The given set of phosphosite references. 
        consurf_d (dict[str, dict[int, float]]): A dictionary mapping ORFs to a dictionary mapping position to ConSurf score.

    Returns:
        list[list[float]]: A list of ConSurf scores for the given references. 
    '''
    result_references = list()
    result_consurf = list()
    for reference in references:
        protein, site = reference.split('_')
        aa = site[0]
        pos = int(site[1:])
        if protein not in consurf_d or pos not in consurf_d[protein]:
            continue
        result_references.append(reference)
        result_consurf.append(consurf_d[protein][pos])
    return result_references, result_consurf

def sample_random_sites(references: set[str], 
                        exclusion: set[str],
                        sequences: dict[str, str], 
                        residues: set[str]) -> set[str]:
    '''
    Given a range of residue types, return residue references sampled within the set of proteins harboring the phosphosites. 
    
    Args: 
        references (set[str]): The given set of phosphosite references. 
        exclusion (set[str]): A set of phosphosite references that have been reported at least once. 
        sequences (dict[str, str]): A dictionary mapping protein name to the sequence. 
        residues (set[str]): Set of residue types that are sampled. 

    Returns:
        set[str]: Residue references sampled. 
    '''
    result_references = []
    # Set of prteins to be sampled
    proteins_unfiltered = set([reference.split('_')[0] for reference in references])
    proteins = set([protein for protein in proteins_unfiltered if protein in sequences])
    for protein in proteins:
        for index, aa in enumerate(sequences[protein]):
            if aa in residues:
                reference = protein + '_' + aa + str(index+1)
                if reference not in exclusion:
                    result_references.append(reference)
    return result_references

def bootstrap(data: list[float], 
              num_iterations=1000, 
              statistic=np.median, 
              seed=None) -> list[float]:
    '''
    Bootstrap resampling to estimate a statistic.

    Args:
        data (list[float]): List of floats. 
        num_interations (int): Number of iterations. 
        statistic (function): Function to apply to each resampled dataset. 
        seed (int): Random seed for reproducibility. 
    '''
    rng = np.random.default_rng(seed)
    data = np.array(data)
    n = len(data)

    bootstraped_statistics = np.array([
        statistic(rng.choice(data, size=n, replace=True)) for _ in range(num_iterations)
    ])

    return bootstraped_statistics

def calculate_ST_component(references: set[str]) -> list[float]:
    '''
    Calculate the fraction of S and T given a set of phosphosite references. 

    Args:
        references (set[str]): The given set of phosphosite references. 

    Returns:
        list[float]: Fraction of S and T, respectively. 
    '''
    S_count = 0
    T_count = 0
    for reference in references: 
        aatype = reference.split('_')[1][0]
        if aatype == 'S':
            S_count += 1
        elif aatype == 'T':
            T_count += 1
        else:
            print(reference)
    print(S_count / (S_count + T_count), T_count / (S_count + T_count))

def sample_proportional(references: set[str], 
                        random_references: set[str], 
                        seed: int) -> set[str]:
    '''
    Sample proteins harboring given references with a probability proportional to the number of phosphosites they have. 
    Use retrieve_references_by_order first to sample in ordered or disordered regions. 

    Args:
        references (set[str]): Given phosphosite references. 
        random_references (set[str]): Random phosphosites to be sampled. 
    '''
    result = []
    if seed is not None: 
        random.seed(seed)
    reference_count = Counter([reference.split('_')[0] for reference in references])
    random_dict = {}
    for random_reference in random_references: 
        protein, site = random_reference.split('_')
        aa = site[0]
        pos = int(site[1:])
        if protein not in random_dict:
            random_dict[protein] = set()
        random_dict[protein].add(random_reference)
    for protein, count in reference_count.items(): 
        if protein in random_dict:
            rs = list(random_dict[protein])
            if len(rs) > count:
                result += random.sample(rs, count)
            else:
                result += rs
    return set(result)

def calculate_consurf_difference_psites(references: set[str],
                                        consurf_d: dict[str, dict[int, float]],
                                        window_size: int = 5) -> list[float]:
    '''
    For a set of residues, calculate difference in ConSurf score between a given residue and adjacent residues (from -2 to +2).

    Args:
        references (set[str]): the set of given residues. 
        consurf_d (dict[str, dict[int, float]]): dictionary mapping protein name to dictionary mapping position to ConSurf score. 
        window_size (int): the size of window for which the difference of ConSurf score is calculated. 

    Returns:
        list[float]: the list of average difference of ConSurf score. 
    '''
    difference_ls = []
    for reference in references:
        diff = calculate_consurf_difference_psite(reference, consurf_d, window_size)
        if diff:
            difference_ls.append(diff)
    return difference_ls

def calculate_consurf_difference_psite(reference: str,  
                                       consurf_d: dict[str, dict[int, float]], 
                                       window_size: int = 5) -> float:
    '''
    Calculate difference in ConSurf score between a given residue and adjacent residues (from -2 to +2).

    Args:
        reference (str): the given residue. 
        consurf_d (dict[str, dict[int, float]]): dictionary mapping protein name to dictionary mapping position to ConSurf score. 
        window_size (int): the size of window for which the difference of ConSurf score is calculated. 

    Returns:
        float: the average difference of ConSurf score. 
    '''
    systematic_name, site = reference.split('_')
    if systematic_name not in consurf_d or systematic_name in ['YKL021C', 'YDR098C', 'YMR112C']:
        return None
    aa = site[0]
    position = int(site[1:])
    if aa not in ['S', 'T']:
        return None
    start = max(position - (window_size - 1) // 2, 1)
    end = min(position + (window_size - 1) // 2, max(consurf_d[systematic_name].keys()))
    if end < start:
        print(systematic_name)
        return None
    consurf_neighbors = []
    for i in range(start, end+1):
        if i == position:
            continue
        consurf_neighbors.append(consurf_d[systematic_name][position] - consurf_d[systematic_name][i])
    return sum(consurf_neighbors) / len(consurf_neighbors)

def retrieve_exposure_references(references: set[str],
                                 resType: str, 
                                 aaInfo: dict[str, dict[int, str]], 
                                 RSAInfo: dict[str, dict[int, str]], 
                                 dRSAInfo: dict[str, dict[int, str]]) -> set[str]:
    '''
    Return references of exposed residues. 

    Args:
        references (set[str]): Set of references for calculation. 
        resType (str): Conditional p-sites, universal p-sites, or random S/T. 
        aaInfo (dict[str, dict[int, str]]): Dictionary recording amino acid type. 
        RSAInfo (dict[str, dict[int, float]]): Dictionary recording residue RSA. 
        dRSAInfo (dict[str, dict[int, float]]): Dictionary recording residue delta RSA. 
    '''
    result_set = set()
    for reference in references: 
        systematic_name, site = reference.split('_')
        if systematic_name in ['YKL021C', 'YDR098C', 'YMR112C']:
            continue
        aa = site[0]
        position = int(site[1:])
        try: 
            if dRSAInfo[systematic_name][position] == 0 and RSAInfo[systematic_name][position] > 0.25:
                result_set.add(reference)
        except KeyError:
            continue
    return result_set

def bootstrap_se(data, n_bootstrap=1000):
    medians = [np.median(np.random.choice(data, size=len(data), replace=True))
               for _ in range(n_bootstrap)]
    return np.std(medians)

def calculate_exposure_consurf(references: set[str], 
                               resType: str, 
                               consurf_d: dict[str, dict[int, float]], 
                               aaInfo: dict[str, dict[int, str]], 
                               RSAInfo: dict[str, dict[int, float]], 
                               dRSAInfo: dict[str, dict[int, float]]) -> list[Union[str, float]]:
    '''
    Calculate ConSurf score by exposure. 

    Args:
        references (set[str]): Set of references for calculation. 
        resType (str): Conditional p-sites, universal p-sites, or random S/T. 
        consurf_d (dict[str, dict[int, float]]): ConSurf score dictionary. 
        aaInfo (dict[str, dict[int, str]]): Dictionary recording amino acid type. 
        RSAInfo (dict[str, dict[int, float]]): Dictionary recording residue RSA. 
        dRSAInfo (dict[str, dict[int, float]]): Dictionary recording residue delta RSA. 
    '''
    interfacial_psite = []
    exposed_psite = []
    buried_psite = []
    for reference in references: 
        systematic_name, site = reference.split('_')
        if systematic_name in ['YKL021C', 'YDR098C', 'YMR112C']:
            continue
        aa = site[0]
        position = int(site[1:])
        try: 
            if dRSAInfo[systematic_name][position] > 0:
                interfacial_psite.append(reference)
            elif RSAInfo[systematic_name][position] > 0.25:
                exposed_psite.append(reference)
            elif RSAInfo[systematic_name][position] <= 0.25:
                buried_psite.append(reference)
        except KeyError:
            continue
    _, consurf_interfacial = retrieve_ConSurf_score(interfacial_psite, consurf_d)
    _, consurf_exposed = retrieve_ConSurf_score(exposed_psite, consurf_d)
    _, consurf_buried = retrieve_ConSurf_score(buried_psite, consurf_d)
    return [('Interfacial', resType, np.median(consurf_interfacial), bootstrap_se(consurf_interfacial)), 
            ('Exposed', resType, np.median(consurf_exposed), bootstrap_se(consurf_exposed)),
            ('Buried', resType, np.median(consurf_buried), bootstrap_se(consurf_buried))]

def calculate_exposure_consurf_alphafold(references, resType, consurf_d, rsa_alphafold, dRSAInfo, window_size=5): 
    interfacial_psite = []
    exposed_psite = []
    buried_psite = []
    for reference in references: 
        systematic_name, site = reference.split('_')
        if systematic_name in ['YKL021C', 'YDR098C', 'YMR112C']:
            continue
        aa = site[0]
        position = int(site[1:])
        try: 
            if dRSAInfo[systematic_name][position] > 0: 
                interfacial_psite.append(reference)
            elif rsa_alphafold[systematic_name][position] > 0.25: 
                exposed_psite.append(reference)
            elif rsa_alphafold[systematic_name][position] <= 0.25: 
                buried_psite.append(reference)
        except KeyError:
            continue
    consurf_interfacial = calculate_consurf_difference_psites(interfacial_psite, consurf_d, window_size)
    consurf_exposed = calculate_consurf_difference_psites(exposed_psite, consurf_d, window_size)
    consurf_buried = calculate_consurf_difference_psites(buried_psite, consurf_d, window_size)
    return [('Interfacial', resType, np.median(consurf_interfacial), bootstrap_se(consurf_interfacial)), 
            ('Exposed', resType, np.median(consurf_exposed), bootstrap_se(consurf_exposed)),
            ('Buried', resType, np.median(consurf_buried), bootstrap_se(consurf_buried))]