#!/usr/bin/env python

import copy
import getopt
import itertools
from itertools import groupby
import logging
import os.path
import statistics as stat
import sys
import time
import random
import re
from collections import Counter

from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import *
import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import (
    silhouette_samples,
    silhouette_score
)
from sklearn.mixture import GaussianMixture
import networkx as nx
import markov_clustering as mcl
from scipy import stats

from .args_processing import process_args
from .helper import (
    read_input_files
)

from .parser import create_parser

def execute(
    tree: str,
    fasta: int,
    support: float,
    occupancy: float
):
    """
    Master execute Function
    -----------------------                                    
    This function executes the main functions and calls other    
    subfunctions to run quest  
    """
    tree, fasta = read_input_files(tree, fasta)
    tree.root_at_midpoint()
    # tree.collapse_all(lambda c: c.confidence is not None and c.confidence < 80)
    taxa = []
    all_tips = []
    for term in tree.get_terminals():
        if term.name.split('|', 1)[0] not in taxa:
            taxa.append(term.name.split('|', 1)[0])
        all_tips.append(term.name)
    
    groups_of_interest = []
    assigned_tips = []
    # skip root, hence [1:]
    for inter in tree.get_nonterminals()[1:]:
        newtree = copy.deepcopy(tree)
        taxa_from_terms = []
        terms = []
        for term in inter.get_terminals():
            taxa_from_terms.append(term.name.split('|', 1)[0])
            terms.append(term.name)
        counts_of_taxa_from_terms = Counter(taxa_from_terms)
        counts = []
        for count in counts_of_taxa_from_terms.values():
            counts.append(count)
        if set([1]) == set(counts) and len(list(set(terms) & set(assigned_tips))) == 0:
            tips_to_prune = [i for i in all_tips + terms if i not in all_tips or i not in terms]
            for tip in tips_to_prune:
                newtree.prune(tip)
            groups_of_interest.append(terms)
            newtree.collapse_all(lambda c: c.confidence is not None and c.confidence < 80)
            # should_get_fasta = True
            if newtree.count_terminals() > (len(taxa)*occupancy):
                Phylo.draw_ascii(newtree)
            for term in inter.get_terminals():
                assigned_tips.append(term.name)            
        elif set([1, 2]) == set(counts) and len(list(set(terms) & set(assigned_tips))) == 0:
            tips_to_prune = [i for i in all_tips + terms if i not in all_tips or i not in terms]
            for tip in tips_to_prune:
                newtree.prune(tip)
            # groups_of_interest.append(terms)
            newtree.collapse_all(lambda c: c.confidence is not None and c.confidence < 80)
            for i in counts_of_taxa_from_terms:
                if counts_of_taxa_from_terms[i] > 1:
                    dups = [e for e in terms if e.startswith(i)]
                    dup_subtrees = []
                    subtree_tips = []
                    for dup in dups:
                        temptree = copy.deepcopy(tree)
                        node_path = temptree.get_path(dup)
                        temp = []
                        for term in node_path[-2].get_terminals():
                            temp.append(term.name)
                        subtree_tips.append(temp)
                    first_set_of_subtree_tips = subtree_tips[0]
                    are_sisters = True
                    for set_of_subtree_tips in subtree_tips[1:]:
                        if first_set_of_subtree_tips != set_of_subtree_tips:
                            are_sisters = False
                        if not are_sisters:
                            break
                    if are_sisters:
                        seq_lengths = dict()
                        for dup in dups:
                            seq_lengths[dup] = len(re.sub('-', '', str(fasta[dup].seq)))
                        longest_seq = max(seq_lengths, key=seq_lengths.get)
                        for seq_name, _ in seq_lengths.items():
                            if seq_name != longest_seq:
                                newtree.prune(seq_name)

                        if newtree.count_terminals() > (len(taxa)*occupancy):
                            Phylo.draw_ascii(newtree)
                        for term in inter.get_terminals():
                            assigned_tips.append(term.name)

def main(argv=None):
    """
    Function that parses and collects arguments              
    """
    # parse and assign arguments
    parser = create_parser()
    args = parser.parse_args()

    # pass to master execute function
    execute(**process_args(args))

if __name__ == "__main__":
    main(sys.argv[1:])
