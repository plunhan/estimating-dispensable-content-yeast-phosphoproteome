'''
Calculate ConSurf scores for S. cerevisiae proteins. 
'''
import os
from lib_rate4site_tools import prepare_rate4site_files, run_rate4site
from pathlib import Path

def main():

    dataDir = Path('../data')

    extDir = dataDir / 'external'

    procDir = dataDir / 'processed'

    rate4siteDir = procDir / 'rate4site_files'

    alignmentDir = procDir / 'MSA'

    ConSurfDir = procDir / 'ConSurf'

    if not rate4siteDir.exists(): 
        os.mkdir(rate4siteDir)

    tree_file = extDir / 'yeast_tree.txt'
    yeast_proteome = extDir / 'Scer.aa'
    rate4sitePath = Path('/Users/plunhan/Desktop/Aim_3/rate4site.3.2.source/sourceMar09/rate4site')

    orfs = prepare_rate4site_files(yeast_proteome, rate4siteDir, alignmentDir)
    run_rate4site(orfs, tree_file, rate4siteDir, ConSurfDir, rate4sitePath)

if __name__ == '__main__':
    main()