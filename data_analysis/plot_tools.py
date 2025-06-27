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

def plot_distribution_consurf(consurf_ls_ls: list[list[float]], 
                              outPath: Union[str, Path], 
                              fmt: str, 
                              label_list: list[str]) -> None: 
    all_data = np.concatenate(consurf_ls_ls)
    x_min, x_max = min(all_data), max(all_data)

    num_bins = 20
    shared_bins = np.linspace(x_min, x_max, num_bins + 1)

    y_max = 0
    for consurf_ls in consurf_ls_ls:
        counts, bins = np.histogram(consurf_ls, bins=20, density=True)
        y_max = max(y_max, max(counts))
    
    fig, axes = plt.subplots(nrows=len(consurf_ls_ls), ncols=1, figsize=(6, 4 * len(consurf_ls_ls)))

    if len(consurf_ls_ls) == 1:
        axes = [axes]

    for i, consurf_ls in enumerate(consurf_ls_ls):
        axes[i].hist(consurf_ls, bins=shared_bins, density=True, edgecolor='black')
        axes[i].set_title(f'{label_list[i]}')
        axes[i].set_xlim(x_min, x_max)
        axes[i].set_ylim(0, y_max * 1.1)
        axes[i].set_ylabel('Probability density')
    plt.xlabel('ConSurf score')
    plt.savefig(outPath, format=fmt, dpi=300)
    plt.close()
'''
def plot_figure2b(ax, data, title, labels):
    sns.boxplot(data=data, ax=ax, showfliers=False)
    ax.set_title(title)
    ax.set_xticklabels(labels, rotation=0)

    y_max = max(max(group) for group in data) + 0.1
    step = 0.05
    count = 0

    for i in range(len(data)):
        for j in range(i+1, len(data)):
            stat, p = ranksums(data[i], data[i+1])
            label = f'p={p:.2e}'
            x1, x2 = i, j
            y = y_max + count * step
            ax.plot([x1, x1, x2, x2], [y, y+0.02, y+0.02, y], lw=1.2, c='k')
            ax.text((x1 + x2) * 0.5, y + 0.025, label, ha='center', va='bottom', fontsize=9)
            count += 1
'''
def melt_data(data, group_labels, ylabel):
    flat_values = sum(data, [])
    group_column = [label for label, group in zip(group_labels, data) for _ in group]
    return pd.DataFrame({'Group': group_column, ylabel: flat_values})

def plot_figure2b(ax, data, title, pairs, labels, ylabel):
    df = melt_data(data, labels, ylabel)
    sns.boxplot(x='Group', y=ylabel, data=df, ax=ax, showfliers=False)
    ax.set_title(title)

    # max_value = df['Value'].max()
    # y_max = max_value + 0.3 * (max_value - df['Value'].min())
    # ax.set_ylim(top=y_max)
    ax.set_ylim(bottom=df[ylabel].min()+0.5)
    
    annotator = Annotator(ax, pairs, data=df, x='Group', y=ylabel)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
    annotator.apply_and_annotate()