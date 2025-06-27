import os
import gzip
import zipfile
import subprocess
from pathlib import Path

def download_initial_data(extDir):
    """
    Get initial input files from SGD, PDB, BioGrid, Ensembl.
    Adapted from Pollet et al. JMB 2022
    Args:
        dir (str): directory to download them into (/data/external directory)
    """
    print ("----------- Download initial data -----------")
    remoteFiles = [
    # ('ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt','pdb_entry_type.txt'),
    # ('ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz','pdb_seqres.txt'),
    # ('ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx','resolu.idx'),
    # ('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz','pdb_chain_uniprot.csv'),                             # Download manually if issues with the REST API (https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
    ('ftp://ftp.ensembl.org/pub/current/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz','Scer.nt'),       # S. cere
    ('ftp://ftp.ensembl.org/pub/current/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz','Scer.aa'),
    ('https://downloads.yeastgenome.org/sequence/fungi/S_uvarum/archive/MIT/orf_dna/orf_genomic.fasta.gz','Sbay.nt'),                           # S. bay
    ('https://downloads.yeastgenome.org/sequence/fungi/S_uvarum/archive/MIT/orf_protein/orf_trans.fasta.gz','Sbay.aa'),
    ('https://downloads.yeastgenome.org/sequence/fungi/S_mikatae/archive/MIT/orf_dna/orf_genomic.fasta.gz','Smik.nt'),                          # S. mik
    ('https://downloads.yeastgenome.org/sequence/fungi/S_mikatae/archive/MIT/orf_protein/archive/orf_trans.20041119.fasta.gz','Smik.aa'),
    ('https://downloads.yeastgenome.org/sequence/fungi/S_paradoxus/archive/MIT/orf_dna/orf_genomic.fasta.gz','Spar.nt'),                        # S. par
    ('https://downloads.yeastgenome.org/sequence/fungi/S_paradoxus/archive/MIT/orf_protein/orf_trans.fasta.gz','Spar.aa'),
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/naumovozyma_castellii_cbs_4309_gca_000237345/cds/Naumovozyma_castellii_cbs_4309_gca_000237345.ASM23734v1.cds.all.fa.gz',"Ncas.nt"),  # N. Cas
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/naumovozyma_castellii_cbs_4309_gca_000237345/pep/Naumovozyma_castellii_cbs_4309_gca_000237345.ASM23734v1.pep.all.fa.gz',"Ncas.aa"),
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/_candida_glabrata_gca_001466525/cds/_candida_glabrata_gca_001466525.ASM146652v1.cds.all.fa.gz',"Cgla.nt"),   #C. Glab
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/_candida_glabrata_gca_001466525/pep/_candida_glabrata_gca_001466525.ASM146652v1.pep.all.fa.gz',"Cgla.aa"),
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/ashbya_gossypii/cds/Ashbya_gossypii.ASM9102v1.cds.all.fa.gz',"Agos.nt"),  #A. Gos
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/ashbya_gossypii/pep/Ashbya_gossypii.ASM9102v1.pep.all.fa.gz',"Agos.aa"),
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/cds/Kluyveromyces_lactis_gca_000002515.ASM251v1.cds.all.fa.gz',"Klac.nt"), #K. Lact
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/pep/Kluyveromyces_lactis_gca_000002515.ASM251v1.pep.all.fa.gz',"Klac.aa"),
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota2_collection/candida_albicans_sc5314_gca_000784635/cds/Candida_albicans_sc5314_gca_000784635.Cand_albi_SC5314_V4.cds.all.fa.gz',"Calb.nt"), # C. Albi
    ('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota2_collection/candida_albicans_sc5314_gca_000784635/pep/Candida_albicans_sc5314_gca_000784635.Cand_albi_SC5314_V4.pep.all.fa.gz',"Calb.aa"),
    ]

    for url, shortName in remoteFiles:
        filePath = extDir / shortName
        if not filePath.is_file():
            cmd_download = ['wget', '-P', str(extDir), url]
            subprocess.run(cmd_download)
            # Rename file to their corresponding shortName, unzip if needed
            outfile = extDir / url.split('/')[-1]
            if (url.split('.')[-1] == "zip"):
                with zipfile.ZipFile(outfile, 'r') as zip_ref:
                    zip_ref.extractall(extDir)
                outfile = Path(str(outfile)[:-4]+".txt")
            if (url.split('.')[-1] == "gz"):
                cmd_gunzip = ['gunzip', str(outfile)]
                subprocess.run(cmd_gunzip)
                outfile = str(outfile)[:-3]
            os.rename(outfile, extDir / shortName)

def main(): 

    dataDir = Path('../data')

    extDir = dataDir / 'external'

    if not dataDir.exists():
        os.mkdir(dataDir)

    if not extDir.exists(): 
        os.mkdir(extDir)

    download_initial_data(extDir)

if __name__ == '__main__':
    main()