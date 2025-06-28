import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn3
from pathlib import Path
from scipy.stats import ranksums
from statannotations.Annotator import Annotator
from typing import Union

def significance_level(p_value: float) -> str:
    if p_value < 0.0001:
        return '****'
    elif p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "n.s."

def plot_count_per_perturbation(count_d: dict[int, set[str]], 
                                outPath: Union[str, Path], 
                                title: str,
                                fmt: str) -> None:
    '''
    Plot count of phosphosites against the number of perturbations the phosphosite is detected. 

    Args:
        count_d (dict[int, set[str]]): A dictionary mapping the number of perturbations to phosphosite references.
        outPath (Union[str, Path]): Path to the output file for the plot. 
        fmt (str): Output file format. 

    Returns:
        None
    '''
    x = sorted(count_d.keys())
    y = [len(count_d[number]) for number in x]
    plt.figure(figsize=(6, 4))
    plt.bar(x, y, color='skyblue')
    plt.xlabel('Number of perturbations in which a phosphosite is detected')
    plt.ylabel('Count of phosphosites')
    plt.title(title)
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_reg_overlap(cond: set[str], 
                     univ: set[str], 
                     reg: set[str],
                     outPath: Union[str, Path], 
                     fmt: str) -> None:
    '''
    Generate a Venn diagram for conditional, universal, and regulated phosphosites. 

    Args:
        cond (set[str]): Set of conditional phosphosites. 
        univ (set[str]): Set of universal phosphosites. 
        reg (set[str]): Set of regulated phosphosites. 

    Returns:
        None
    '''
    venn3([cond, univ, reg], 
          set_labels=('Conditional p-sites', 'Universal p-sites', 'Regulated p-sites'), 
          set_colors=('orange', 'blue', 'red'), alpha=0.7)
    plt.title("Overlap between conditional, universal, and regulated p-sites")
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_comparison_ConSurf(consurf_1: list[float],
                            consurf_2: list[float],
                            label_1: str,
                            label_2: str,
                            outPath: Union[str, Path],
                            fmt: str) -> None:
    '''
    Plot comparison between two groups of ConSurf scores. 

    Args:
        consurf_1 (list[float]): ConSurf scores of the first group. 
        consurf_2 (list[float]): ConSurf scores of the second group.
        label_1 (str): Label of the first group. 
        label_2 (str): Label fo the second group. 
        outPath (Union[str, Path]): Path to the output file for the plot. 
        fmt (str): Output file format. 

    Returns:
        None
    '''
    data = consurf_1 + consurf_2
    labels = [label_1] * len(consurf_1) + [label_2] * len(consurf_2)

    plt.figure(figsize=(6, 5))
    ax = sns.boxplot(x=labels, y=data, showfliers=False)
    plt.xlabel('Group of phosphosites')
    plt.ylabel('ConSurf score')
    y_max = max(data) + 0.1
    x1, x2 = 0, 1
    plt.plot([x1, x1, x2, x2], 
             [y_max, y_max + 0.05, y_max + 0.05, y_max],
             lw=1.5, 
             c='black')
    stat, p_value = ranksums(consurf_1, consurf_2)
    sig = significance_level(p_value)
    plt.text((x1 + x2) / 2,
             y_max + 0.07,
             sig,
             ha='center', 
             va='bottom',
             fontsize=12)
    plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def bootstrap_se(data, n_bootstrap=1000):
    medians = [np.median(np.random.choice(data, size=len(data), replace=True))
               for _ in range(n_bootstrap)]
    return np.std(medians)

def plot_consurf_distribution(data_lists: list[list[float]], 
                              group_names: list[str], 
                              outPath: Union[str, Path], 
                              fmt: str, 
                              urge_positive=False) -> None:

    medians = [np.median(d) for d in data_lists]
    ses = [bootstrap_se(d) for d in data_lists]

    fig, ax = plt.subplots(figsize=(6, 5))
    x_pos = np.arange(len(data_lists))

    ax.bar(x_pos, medians, yerr=ses, capsize=10, width=0.5, edgecolor='black', linewidth=1.2)

    pairs = [(0, 1), (1, 2), (0, 2)]
    if urge_positive:
        height_buffer = 0.1
    else:
        height_buffer = max(medians) + max(ses) * 1.5
    step = max(ses) * 2
    for i, (i1, i2) in enumerate(pairs):
        stat, pvalue = ranksums(data_lists[i1], data_lists[i2])
        sig = significance_level(pvalue)
        y = height_buffer
        if i == 2: 
            y = height_buffer + step
        x1, x2 = x_pos[i1], x_pos[i2]
        ax.plot([x1+0.02, x1+0.02, x2-0.02, x2-0.02], [y, y + 0.02, y + 0.02, y], lw=1.2, c='black')
        ax.text((x1 + x2) / 2, y + 0.025, sig, ha='center', va='bottom', fontsize=10)

    if urge_positive: 
        ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_names)
    ax.set_ylabel("Median ConSurf score")
    ax.set_title("Comparison of median ConSurf score")
    # plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_consurf_difference(data_lists: list[list[float]], 
                            group_names: list[str], 
                            outPath: Union[str, Path], 
                            fmt: str) -> None:

    medians = [-np.median(d) for d in data_lists]
    ses = [bootstrap_se(d) for d in data_lists]

    fig, ax = plt.subplots(figsize=(6, 5))
    x_pos = np.arange(len(data_lists))

    ax.bar(x_pos, medians, yerr=ses, capsize=10, width=0.5, edgecolor='black', linewidth=1.2)

    pairs = [(0, 1), (1, 2), (0, 2)]
    height_buffer = max(medians) + max(ses) * 1.5
    step = max(ses) * 2
    for i, (i1, i2) in enumerate(pairs):
        stat, pvalue = ranksums(data_lists[i1], data_lists[i2])
        sig = significance_level(pvalue)
        y = height_buffer
        if i == 2: 
            y = height_buffer + step
        x1, x2 = x_pos[i1], x_pos[i2]
        ax.plot([x1+0.05, x1+0.05, x2-0.05, x2-0.05], [y, y + 0.02, y + 0.02, y], lw=1.2, c='black')
        ax.text((x1 + x2) / 2, y + 0.025, sig, ha='center', va='bottom', fontsize=10)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_names)
    ax.set_ylabel("Median of difference in ConSurf score between\nadjacent residues and p-site")
    ax.set_title("Difference in ConSurf score between adjacent residues and p-site")
    # plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_consurf_distribution_against_perturbations(data_lists: list[list[str]], 
                                                    group_names: list[str], 
                                                    outPath: Union[str, Path], 
                                                    fmt: str) -> None:
    '''
    Plot distribution of ConSurf score against the number of perturbations. 
    Perturbation bins: 1, 2-4, 5-7, 8-10

    Args: 
        references (set[str]): Set of references. 
        consurf (dict[str, dict[int, float]]): ConSurf dictionary. 
        outPath (Union[str, Path]): Path to output file. 
        fmt (str): Format of the output file. 
    
    Returns:
        None
    '''
    medians = [np.median(d) for d in data_lists]
    ses = [bootstrap_se(d) for d in data_lists]

    fig, ax = plt.subplots(figsize=(6, 5))
    x_pos = np.arange(len(data_lists))

    ax.bar(x_pos, medians, yerr=ses, capsize=10, width=0.5, edgecolor='black', linewidth=1.2)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_names)
    ax.set_ylabel("Median ConSurf score")
    ax.set_title("Distribution of ConSurf score against number of perturbations")
    # plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()