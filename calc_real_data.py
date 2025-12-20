from Bio import Phylo
import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
from synteny_tools import findSyntenyReal1
from synteny_tools import findSyntenyReal2


def get_real_genomes_from_cc(dataset):
    df = pd.read_csv(f"ATGC{dataset}/atgc.cc.csv", names=[
        "gene_ID", "genome_ID", "protein_ID", "protein_length",
        "atgc_cog_footprint", "atgc_cog_footprint_length",
        "cog_ID", "cls_ID", "match_class"
    ])

    df = df.dropna(subset=["cls_ID"])

    return df.groupby("genome_ID")["cls_ID"].apply(list).to_dict()

def convert_to_numeric(cls_list):
    return [int(cls.replace("cls.", "")) for cls in cls_list]

tree = Phylo.read("ATGC0070/atgc.iq.r.tre", "newick")
genome_to_genes = get_real_genomes_from_cc("0070")
print(genome_to_genes)
leaf_names = [leaf.name for leaf in tree.get_terminals()]

# Dict to store block lengths per genome pair

def get_pair_blocks():
    real_pairwise_blocks = {}
    for genome1, genome2 in combinations(leaf_names, 2):
        if genome1 not in genome_to_genes or genome2 not in genome_to_genes:
            continue  # skip missing genomes

        g1 = convert_to_numeric(genome_to_genes[genome1])
        g2 = convert_to_numeric(genome_to_genes[genome2])

        blocks = findSyntenyReal2(g1.copy(), g2.copy())
        block_lengths = [b[4] for b in blocks]
        real_pairwise_blocks[(genome1, genome2)] = block_lengths
    return real_pairwise_blocks

def get_pair_blocks2():
    real_pairwise_blocks = {}
    for genome1, genome2 in combinations(leaf_names, 2):
        if genome1 not in genome_to_genes or genome2 not in genome_to_genes:
            continue  # skip missing genomes

        g1 = convert_to_numeric(genome_to_genes[genome1])
        g2 = convert_to_numeric(genome_to_genes[genome2])

        blocks = findSyntenyReal1(g1.copy(), g2.copy())
        block_lengths = [b[4] for b in blocks]
        real_pairwise_blocks[(genome1, genome2)] = block_lengths
    return real_pairwise_blocks

print(get_pair_blocks2() == get_pair_blocks())
# #optional plot
# real_blocks = get_pair_blocks()

# all_lengths = [length for blocklist in real_blocks.values() for length in blocklist]
# plt.figure(figsize=(10, 5))
# plt.hist(all_lengths, bins=range(1, max(all_lengths) + 2), color="skyblue", edgecolor="black")
# plt.title("Distribution of Synteny Block Lengths (Real Genome Pairs)")
# plt.xlabel("Synteny Block Length")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.show()

# data_setup_from_pty.py
# calc_real_blocks.py


