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
                              title: str,
                              step: float,
                              ylims: tuple[float],
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
    ax.set_ylim(ylims)
    ax.set_title(f"Comparison of median ConSurf score in {title} regions")
    # plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_consurf_distribution_separate(data_lists: list[list[float]], 
                                       group_names: list[str], 
                                       outPath: Union[str, Path], 
                                       fmt: str, 
                                       title: str,
                                       residue: str,
                                       urge_positive=False) -> None:
    medians = [np.median(d) for d in data_lists]
    ses = [bootstrap_se(d) for d in data_lists]

    fig, ax = plt.subplots(figsize=(6, 3))
    x_pos = np.arange(len(data_lists))

    ax.bar(x_pos, medians, yerr=ses, capsize=10, width=0.5, edgecolor='black', linewidth=1.2)

    pairs = [(0, 1), (1, 2), (0, 2)]
    if urge_positive:
        height_buffer = 0.1
    else:
        height_buffer = max(medians) + max(ses) * 1.5

    if residue == 'S':
        step = max(ses) * 1.5
        ax.set_ylim(0, 0.5)
    elif residue == 'T':
        step = 0.1
        ax.set_ylim(0, 1.0)
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
    ax.set_title(title)
    # plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_consurf_difference(data_lists: list[list[float]], 
                            group_names: list[str], 
                            outPath: Union[str, Path], 
                            ylims: tuple[float],
                            fmt: str) -> None:

    medians = [np.median(d) for d in data_lists]
    ses = [bootstrap_se(d) for d in data_lists]

    fig, ax = plt.subplots(figsize=(6, 5))
    x_pos = np.arange(len(data_lists))

    ax.bar(x_pos, medians, yerr=ses, capsize=10, width=0.5, edgecolor='black', linewidth=1.2)

    pairs = [(0, 1), (1, 2), (0, 2)]
    height_buffer = 0.025
    step = max(ses) * 2
    for i, (i1, i2) in enumerate(pairs):
        stat, pvalue = ranksums(data_lists[i1], data_lists[i2])
        sig = significance_level(pvalue)
        y = height_buffer
        if i == 2: 
            y = height_buffer + step
        x1, x2 = x_pos[i1], x_pos[i2]
        ax.plot([x1+0.025, x1+0.025, x2-0.025, x2-0.025], [y, y + 0.02, y + 0.02, y], lw=1.2, c='black')
        ax.text((x1 + x2) / 2, y + 0.025, sig, ha='center', va='bottom', fontsize=10)

    ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_names)
    ax.set_ylabel("Median ConSurf score of p-site relative to ajacent residues")
    ax.set_title("ConSurf score of p-site relative to adjacent residues")
    ax.set_ylim(ylims)
    # plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_consurf_distribution_against_perturbations(data_lists: list[list[str]], 
                                                    group_names: list[str], 
                                                    outPath: Union[str, Path], 
                                                    fmt: str, 
                                                    ylims: tuple[float], 
                                                    title) -> None:
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
    ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_names)
    ax.set_xlabel("Number of perturbations")
    ax.set_ylabel("Median ConSurf score")
    ax.set_ylim(ylims)
    ax.set_title(title)
    # plt.tight_layout()
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()

def plot_consurf_exposure(df: pd.DataFrame, 
                          outPath: Union[str, Path], 
                          figFmt: str) -> None:
    plt.figure(figsize=(7, 5))
    ax = sns.barplot(
        data=df, 
        x='Exposure', 
        y='Median',
        hue='Type', 
        palette='Set2',
        ci=None
    )

    for bar in ax.patches:
        bar.set_edgecolor('black')
        bar.set_linewidth(1.2)

    for bar, (_, row) in zip(ax.patches, df.iterrows()):
        x = bar.get_x() + bar.get_width() / 2
        y = bar.get_height()
        ax.errorbar(x, y, yerr=row['Standard error'], fmt='none', 
                    ecolor='black', capsize=5, elinewidth=1)

    ax.axhline(0, color='black', linestyle='-', linewidth=1)

    ax.set_title("Distribution of evolutionary conservation for residues with different exposure")
    ax.set_xlabel("Residue exposure")
    ax.set_ylabel("Median ConSurf score")
    plt.legend(title="Residue type")
    plt.savefig(outPath, dpi=300, format=figFmt)
    plt.close()