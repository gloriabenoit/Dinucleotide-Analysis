"""Module to analyze a fasta file's dinucleotides.

Usage:
======
    python GloriaBENOIT_DinucleotideAnalysis.py argument1 ... argumentN

    argument1 through N : caracter chain for a zipped fasta file
"""

__authors__ = "Gloria Benoit"
__date__ = "2023-04-23"

import math
import sys
import os
import gzip

import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
import scipy.stats as stats
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage

DINUCL = ["AA", "AC", "AG", "AT", "CA", "CC", "CG",
          "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
BASES = ['A', 'T', 'C', 'G']
SPECIES = ["S. pneumoniae", "E. coli", "A. thaliana", "S. cerevisiae",
           "P. dumerilii", "C. elegans", "D. melanogaster", "A. rubens",
           "P. marinus", "S. bombifrons", "P. cristatus", "M. musculus"]


def initialize_dictio(liste):
    """Initialize a dictionary based on list values.

    This function creates a dictionary whose keys are the
    list's values initialized at 0.

    Parameters
    ----------
    liste : list
        The list of the values we want to store.

    Returns
    -------
    dict
        A new dictionary.
    """
    dictio = {}
    for items in liste:
        dictio[items] = 0
    return dictio


def show_dictio(dictio):
    """Show the content of the dictionary.

    Parameters
    ----------
    dictio : dict
        The dictionary you want to show.
    """
    for key, val in dictio.items():
        print(f"{key} : {val}")


def add_count(value, dictio):
    """Count the value if it's already in the dictionary.

    Parameters
    ----------
    value : any
        The value you want to count.
    dictio : dict
        The dictionary to store the counts.

    Returns
    -------
    dict
        A dictionary where the value is incremented by one.
    """
    if value in dictio:
        dictio[value] += 1
    return dictio


def dinucleotide_count(nb_dinucl, sequence, reverse=True):
    """Count each dinucleotide in a sequence.

    This function counts each dinucleotide on both strand.

    Parameters
    ----------
    nb_dinucl : dict
        A dictionary whose keys are each dinucleotide.
    sequence : str
        The DNA sequence.
    reverse : bool, optional
        If True (default), takes into account both DNA strand.
        If False, takes only the first strand.

    Returns
    -------
    dict
        A dictionary containing all occurrences of each
        dinucleotide found in the sequence.
    """
    for i in range(len(sequence) - 1):
        coding = Seq(sequence[i] + sequence[i + 1])
        add_count(coding, nb_dinucl)
        if (reverse):
            add_count(coding.reverse_complement(), nb_dinucl)
    return nb_dinucl


def expected_base_count(expected_base, sequence):
    """Count each base in a sequence.

    Parameters
    ----------
    expected_base : dict
        A dictionary whose keys are each base.
    sequence : str
        The DNA sequence.

    Returns
    -------
    dict :
        A dictionary containing all occurrences of each
        base found in the sequence.
    """
    for base in BASES:
        nb = sequence.count(base)
        expected_base[base] += nb
        expected_base[(Seq(base)).reverse_complement()] += nb
    return expected_base


def expected_dinucleotide_frequency(base_freq):
    """Calculate the expected frequency of each dinucleotide.

    Parameters
    ----------
    base_freq : dict
        A dictionary whose keys are each base and values
        are their occurrences found in a sequence.

    Returns
    -------
    dict :
        A dictionary containing all expected frequencies
        of each dinucleotide found in the sequence.
    """
    expected_nb_dinucl = initialize_dictio(DINUCL)
    for key in expected_nb_dinucl:
        frequency = base_freq[key[0]] * base_freq[key[1]]
        expected_nb_dinucl[key] = frequency
    return expected_nb_dinucl


def freq_to_count(dictio, length):
    """Transform occurrences in frequencies.

    Parameters
    ----------
    dictio : dict
        A dictionary whose values are the occurrences.
    length : int
        The value to divide the occurrences by.

    Returns
    -------
    dict
        A dictionary whose values are the frequencies.
    """
    new_dictio = {}
    for key, val in dictio.items():
        new_dictio[key] = val / length
    return new_dictio


def bootstrap_chi2(obs, exp, nb_val=1000, nb_loop=1000):
    """Average of calculated chi2 using bootstrap.

    Parameters
    ----------
    obs : float
        The observed CG frequency.
    exp : float
        The expected CG frequency.
    nb_val : int, optional
        The sample size (default = 1000).
    nb_loop : int, optional
        The number of bootstrap performed (default = 1000).

    Returns
    -------
    bool
        The significance of the difference between CG frequencies.
    """
    sumChi2 = []
    d_theo = int(exp * nb_val)
    nb_signif = 0

    val = np.array([1, 0])  # 1 is CG and 0 is other
    prob = np.array([obs, 1-obs])

    for i in range(nb_loop):
        ech = np.random.choice(val, size=nb_val, replace=True, p=prob)
        d_obs = np.count_nonzero(ech == 1)
        table = [[d_obs, nb_val - d_obs],
                 [d_theo, nb_val - d_theo]]
        chi_squared_result = stats.chisquare(table[0], table[1])
        sumChi2.append(chi_squared_result[0])
        p = chi_squared_result[1]
        if p < 0.05:
            nb_signif += 1

    # Mean of the chi-squared values
    mean_chi2 = sum(sumChi2) / nb_loop

    # Degrees of freedom
    ddl = nb_val - 1

    # Global p-value of all tests
    p_val = stats.distributions.chi2.sf(mean_chi2, ddl)

    if p < 0.05:
        return True
    else:
        return False


def save_compar_frequency(observed, expected, organism, signif, db=False):
    """Save a bar plot comparing two frequencies.

    Parameters
    ----------
    observed : dict
        A dictionary whose keys are each dinucleotide, and values are
        their observed frequencies.
    expected : dict
        A dictionary whose keys are each dinucleotide, and values are
        their expected frequencies.
    organism : str
        The name of the organism.
    signif : bool
        Indicates whether the CG frequency is significatively different
        between groups.
    db : bool, optional
        If False (default), compares observed and expected frequencies.
        If True, compares observed frequencies with double strand
        and single strand.

    Returns
    -------
    fig
        The figure obtained.
    """
    # Details
    width = 0.35
    colors = ["#cdb4db", "#ffafcc"]
    labels = ["observed", "expected"]
    if (db):
        colors = ["#f4a261", "#e76f51"]
        labels = ["Both strand", "Single strand"]

    # Barplot values
    fig, ax = plt.subplots()
    x = np.arange(len(observed))
    obs = ax.bar(x - width/2, observed.values(), width,
                 color=colors[0], label=labels[0])
    exp = ax.bar(x + width/2, expected.values(), width,
                 color=colors[1], label=labels[1])

    # Barplot details
    ax.set_xticks(x, observed.keys())
    ax.set(ylim=(0, 0.14))
    ax.set_xlabel("Dinucleotide")
    ax.set_ylabel("Frequency")
    ax.set_title(f"Dinucleotide frequency of {organism}")
    ax.legend()

    if (signif):
        # Star above CG frequency
        ax.scatter(6, max(observed["CG"], expected["CG"])+0.005,
                   color="black", marker='*', s=25)

    if (db):
        # Comparison between both strand and single strand
        plt.savefig(f"../results/{organism}_db_sb_freq_dinucl.png", dpi=300,
                    bbox_inches="tight")

    if (not db):
        # Comparison between observed and expected frequencies
        plt.savefig(f"../results/{organism}_obs_exp_freq_dinucl.png", dpi=300,
                    bbox_inches="tight")
        return fig


def save_clustering(names, CG_values):
    """Save a dendrogram of species based on their CG frequencies.

    Parameters
    ----------
    names : list
        A list of all species.
    CG_values : list
        A list of all CG frequencies.
    """
    # Distance matrix
    dists = pdist(np.array(CG_values).reshape(-1, 1))

    # Clustering
    Z = linkage(dists, "ward")

    # Dendrogram values
    fig, ax = plt.subplots(figsize=(10, 5))
    dendrogram(Z, labels=names, ax=ax, orientation="left")

    # Dendrogram details
    plt.title("Species clustering based on CG frequency")
    plt.xlabel("Species")
    plt.ylabel("Distances")
    plt.savefig(f"../results/species_clustering.png", dpi=300,
                bbox_inches="tight")


def save_multiple_fig(fig_list):
    """Save multiple figures displayed at once.

    Parameters
    ----------
    fig_list : list
        A list of matplotlib figures.
    """
    # Space details
    if (len(fig_list) <= 3):
        ncols = 2
        nrow = 1
    ncols = math.ceil(len(fig_list) / 3)
    if (ncols == 1):
        ncols += 1
    nrows = 3
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(ncols*4, nrows*3))

    # Plotting all figs
    for i, ax_row in enumerate(axes):
        for j, ax in enumerate(ax_row):
            fig_index = i * ncols + j
            if fig_index < len(fig_list):
                ax.imshow(fig_list[fig_index].canvas.renderer.buffer_rgba(),
                          aspect="auto")
                ax.axis("off")
            else:
                # Hide the subplot if there's no corresponding figure
                ax.set_visible(False)

    # Complete plot details
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    plt.savefig(f"../results/all_figures.png", dpi=300,
                bbox_inches="tight")


def save_results(name, seq, base, unresolved, obs, exp, signi):
    """Save results' information.

    Parameters
    ----------
    name : str
        The name of the species.
    seq : int
        The number of sequences in the genome.
    base : int
        The number of bases in the genome.
    unresolved : float
        The percentage of unresolved bases.
    obs : dict
        A dictionary whose keys are each dinucleotide, and values are
        their observed frequencies.
    exp : dict
        A dictionary whose keys are each dinucleotide, and values are
        their expected frequencies.
    signif : bool
        Indicates whether the CG frequency is significatively different
        between groups.
    """
    # Create/ write into result file
    with open("results.txt", "a") as filout:
        if (unresolved < 0.01):
            filout.write(f"{name} - {seq} sequences, {base} bases, "
                         f"<0.01% unresolved\n")
        else:
            filout.write(f"{name} - {seq} sequences, {base} bases, "
                         f"{unresolved:.2f}% unresolved\n")
        filout.write("Observed frequencies of each dinucleotide :\n")
        for key, val in obs.items():
            filout.write(f"{key} : {val}\n")
        filout.write("Expected frequencies of each dinucleotide :\n")
        for key, val in exp.items():
            filout.write(f"{key} : {val}\n")
        if (signi):
            filout.write("Observed and expected CG frequencies "
                         "are significatively different (T).\n")
        else:
            filout.write("Observed and expected CG frequencies "
                         "are not significatively different (F).\n")
        filout.write("\n")


def get_dinucleotide_analysis(list_names, list_obs, list_exp,
                              list_signif, filename="results.txt"):
    """Retrieve dinucleotide analysis from result file.

    Parameters
    ----------
    list_names : list of str
        A list of names to add to.
    list_obs : list of dict
        A list of observed dinucleotide frequencies to add to.
    list_exp : list of dict
        A list of expected dinucleotide frequencies to add to.
    list_signif : list of bool
        A list of significance to add to.
    filename: str
        Name of the result file (default = "results.txt")

    Returns
    -------
    list of str
        Updated list of names.
    list of dict
        Updated list of observed dinucleotide frequencies.
    list of dict
        Updated list of expected dinucleotide frequencies.
    list of bool
        updated list of significance.
    """
    flag = -1  # To specify if frequency is observed or expected
    dictio_obs = {}
    dictio_exp = {}
    with open(filename, "r") as result:
        for ligne in result:
            if (ligne.endswith("unresolved\n")):
                name = ligne.split("-")[0]
                name = name.strip()
                list_names.append(name)
            elif (ligne.startswith("Observed")):
                if (ligne.endswith(":\n")):
                    # Next lines will be observed
                    flag += 1
                elif (ligne.endswith("(F).\n")):
                    list_signif.append(False)
                elif (ligne.endswith("(T).\n")):
                    list_signif.append(True)
            elif (ligne.startswith("Expected")):
                # Next lines will be expected
                flag = -1
            elif (ligne == "\n"):
                list_obs.append(dictio_obs)
                list_exp.append(dictio_exp)
                dictio_obs = {}
                dictio_exp = {}
            else:
                if (flag == 0):
                    key = ligne.split(" : ")[0]
                    val = float(ligne.split(" : ")[1])
                    dictio_obs[key] = val
                else:
                    key = ligne.split(" : ")[0]
                    val = float(ligne.split(" : ")[1])
                    dictio_exp[key] = val
    return list_names, list_obs, list_exp, list_signif


def new_dinucleotide_analysis(list_names, list_obs, list_exp,
                              list_signif, filenames, names,
                              strand=False, verbose=False):
    """Analyze dinucleotides in a zipped fasta file.

    Parameters
    ----------
    list_names : list of str
        A list of names to add to.
    list_obs : list of dict
        A list of observed dinucleotide frequencies to add to.
    list_exp : list of dict
        A list of expected dinucleotide frequencies to add to.
    list_signif : list of bool
        A list of significance to add to.
    filenames : list
        A list of all files.
    names : list
        A list of all species in file.
    strand : bool
        If True, analysis of single and double strand frequencies
        conducted on first species.
    verbose : bool, optional
        If True, explains in details what is happening.

    Returns
    -------
    list of str
        Updated list of names.
    list of dict
        Updated list of observed dinucleotide frequencies.
    list of dict
        Updated list of expected dinucleotide frequencies.
    list of bool
        Updated list of significance.
    """
    # Parse each sequence of the genome
    for i in range(len(filenames)):
        if not os.path.exists(filenames[i]):
            continue
            print(f"{filenames[i]} does not exist in the working "
                  "directory, it has not been analyzed.")

        # Number of sequences and total length of the genome
        cpt, length = 0, 0
        dinucl_obs_count = initialize_dictio(DINUCL)
        expected_base = initialize_dictio(BASES)

        # Analysis of single strand for first species
        if (strand):
            if (verbose):
                print("We will compare double and single strand "
                      f"dinucleotide frequencies for {names[i]}.")
            if (i == 0):
                noreverse = initialize_dictio(DINUCL)

        if (verbose):
            print(f"\n{names[i]} :")
        with gzip.open(filenames[i], "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                cpt += 1
                if ((cpt % 10) == 0):
                    if (verbose):
                        print(f"Sequence {cpt}")
                sequence = str(record.seq)
                length += len(sequence)

                # Dinucleotide and base counts
                dinucl_obs_count = dinucleotide_count(dinucl_obs_count,
                                                      sequence)
                expected_base = expected_base_count(expected_base,
                                                    sequence)

                # Analysis of single strand for first species
                if (strand):
                    if (i == 0):
                        noreverse = dinucleotide_count(noreverse,
                                                       sequence,
                                                       reverse=False)

        if (verbose):
            print(f"We counted {cpt} sequences for {names[i]}. "
                  f"The total length of its genome is {length} bases.")

        # Calculate the percentage of unknown bases in the genome
        len_known = sum(expected_base.values())/2
        freq_unknown = ((length - len_known) / length) * 100
        if (verbose):
            if (freq_unknown < 0.01):
                print(f"There is less than 0.01% of unsolved bases.")
            else:
                print(f"There is {freq_unknown:.2f}% of unsolved bases.")

        # Obtain frequencies from occurrences
        dinucl_obs_freq = freq_to_count(dinucl_obs_count, (2*len_known-2))
        expected_base = freq_to_count(expected_base, (2*len_known))
        dinucl_exp_freq = expected_dinucleotide_frequency(expected_base)

        # Analysis of single strand for first species
        if (strand):
            if (i == 0):
                noreverse = freq_to_count(noreverse, (len_known-1))
                save_compar_frequency(dinucl_obs_freq, noreverse,
                                      names[i], signif=False, db=True)
                if (verbose):
                    print("We have saved the comparison of "
                          "dinucleotide frequencies between double "
                          f" and single strand for {names[i]}.")

        if (verbose):
            print("We have obtained the observed and expected "
                  f"frequencies for {names[i]}.\n"
                  "Observed frequencies of each dinucleotide :")
            show_dictio(dinucl_obs_freq)
            print("Expected frequencies of each dinucleotide :")
            show_dictio(dinucl_exp_freq)

        # Statistic test for homogeneity of CG counts
        signif = bootstrap_chi2(dinucl_obs_freq["CG"], dinucl_exp_freq["CG"])
        if (verbose):
            if (signif):
                print("Observed and expected CG frequencies are "
                      "significatively different.")
            else:
                print("Observed and expected CG frequencies are "
                      "not significatively different.")

        # Saving results
        if (verbose):
            print("The results are saved in 'results.txt'.")
        save_results(names[i], cpt, length, freq_unknown,
                     dinucl_obs_freq, dinucl_exp_freq, signif)

        # Adding all values to lists
        list_names.append(names[i])
        list_obs.append(dinucl_obs_freq)
        list_exp.append(dinucl_exp_freq)
        list_signif.append(signif)

    return list_names, list_obs, list_exp, list_signif


def full_dinucleotide_analysis(list_names, list_obs, list_exp,
                               list_signif):
    """Full dinucleotide analysis from lists.

    Parameters
    ----------
    list_names : list of str
        A list of names to add to.
    list_obs : list of dict
        A list of observed dinucleotide frequencies to add to.
    list_exp : list of dict
        A list of expected dinucleotide frequencies to add to.
    list_signif : list of bool
        A list of significance to add to.
    filename: str
        Name of the result file (default = "results.txt")
    """
    figures = []
    CG_values = []
    for i in range(len(list_names)):
        # Figure of observed and expected frequencies
        fig = save_compar_frequency(list_obs[i], list_exp[i],
                                    list_names[i], list_signif[i])
        figures.append(fig)
        if (verbose):
            print("We have saved the comparison of observed and expected"
                  f" dinucleotide frequencies for {list_names[i]}.")

        # Saving total CG frequency
        CG_values.append(list_obs[i]["CG"])

    if (len(figures)) > 1:
        # Hierarchal clustering based on CG frequencies
        save_clustering(list_names, CG_values)
        if (verbose):
            print("We have saved the hierarchal clustering of species based"
                  " on their CG frequencies.")

        # All results shown at once
        save_multiple_fig(figures)
        if (verbose):
            print("We have saved the comparison of observed and expected"
                  " dinucleotide frequencies for all species.")


if __name__ == "__main__":
    # Check for values
    if (len(sys.argv) == 1) and not (os.path.exists("results.txt")):
        sys.exit("ERROR : please enter at least one fasta file "
                 "to start the analysis.")

    print("All figures will be saved in the working directory.")
    print("Please make sure you changed 'SPECIES' according to"
          " the files you'll be analyzing.")
    verbose = 0
    while ((verbose != 'Y') & (verbose != 'N')):
        verbose = input("Do you want to activate verbose mode?"
                        "\nY/N\n")
    if (verbose == 'Y'):
        verbose = True
    else:
        verbose = False

    list_names = []
    list_obs = []
    list_exp = []
    list_signif = []

    # Check for previous results
    if os.path.exists("results.txt"):
        res = get_dinucleotide_analysis(list_names, list_obs,
                                        list_exp, list_signif)
        list_names = res[0]
        list_obs = res[1]
        list_exp = res[2]
        list_signif = res[3]
        if (verbose):
            print("Previous results were recovered.")

    # Get new results
    if (len(sys.argv) > 1):
        files = sys.argv[1:]
        names = SPECIES[:len(files)]
        print("Analysis of new input species.")
        res = new_dinucleotide_analysis(list_names, list_obs,
                                        list_exp, list_signif,
                                        files, names,
                                        verbose=verbose)
        list_names = res[0]
        list_obs = res[1]
        list_exp = res[2]
        list_signif = res[3]

    # Complete analysis
    full_dinucleotide_analysis(list_names, list_obs, list_exp,
                               list_signif)
    if (verbose):
        print("The analysis is over.")