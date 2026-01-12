from Bio import Phylo
from io import StringIO
from ete3 import Tree
import random, numpy as np
import pandas as pd

# function to traverse tree, might be useful later
def traverse(clade, parent_name=""):
    name = clade.name or parent_name  
    print(f"Node {name} (branch length from parent = {clade.branch_length})")
    for child in clade.clades:
        traverse(child, parent_name=name)

def evolve_genome_branch(genome, branch_length, gain_rate, loss_rate, inv_rate):
    genome = genome.copy()  
    # total num of events
    # total_rate = gain_rate + loss_rate + inv_rate
    # N = np.random.poisson(total_rate * branch_length) 
    
    # for _ in range(N):
    #     r = random.random() * total_rate
    #     if r < loss_rate:
    #         # gene loss
    #         if genome:
    #             idx = random.randrange(len(genome))
    #             genome.pop(idx)
    #     elif r < loss_rate + gain_rate:
    #         # gene gain
    #         new_gene = f"newGene{random.randrange(10000)}"
    #         idx = random.randrange(len(genome)+1)
    #         genome.insert(idx, new_gene)
    #     else:
    #         # inversion
    #         if len(genome) >= 2:
    #             i = random.randrange(len(genome))
    #             j = random.randrange(i, len(genome))
    #             segment = genome[i:j+1]
    #             segment.reverse()
    #             genome[i:j+1] = segment

    # changed code to now handle each event independently; need to confirm this is acurate to delete old code
    gains = np.random.poisson(gain_rate * branch_length)
    losses = np.random.poisson(loss_rate * branch_length) 
    inversions = np.random.poisson(inv_rate * branch_length)
    for _ in range(gains):
        new_gene = f"newGene{random.randrange(10000)}"
        idx = random.randrange(len(genome)+1)
        genome.insert(idx, new_gene)
    for _ in range(losses):
        if genome:
            idx = random.randrange(len(genome))
            genome.pop(idx)
    for _ in range(inversions):
        if len(genome) >= 2:
            i = random.randrange(len(genome))
            j = random.randrange(i, len(genome))
            segment = genome[i:j+1]
            segment.reverse()
            genome[i:j+1] = segment
    return genome

def evolve_genome(tree, root_genome, per_gene_gain_rate, per_gene_loss_rate, per_gene_inv_rate):
    genomes = {}
    num_genes = len(root_genome)

    median_path_length = median_root_to_leaf_lengths(tree)
    branch_gain_rate = (num_genes * per_gene_gain_rate) / median_path_length
    branch_loss_rate = (num_genes * per_gene_loss_rate) / median_path_length
    branch_inv_rate = (num_genes * per_gene_inv_rate) / median_path_length

    genomes[tree.root] = root_genome.copy() #now using same root as real data

    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            parent_genome = genomes[parent]
            child_genome = evolve_genome_branch(parent_genome, child.branch_length,
                                                branch_gain_rate, branch_loss_rate, branch_inv_rate)
            genomes[child] = child_genome
    return genomes

def median_root_to_leaf_lengths(tree):
    return np.median([tree.distance(tree.root, leaf) for leaf in tree.get_terminals()])


#FIX THIS to account for inversions ie gene1, gene2, gene3 and gene3 gene2 and gene1 are synteny blocks; inversion do not affect
########
def synteny_blocks(genome1, genome2):

    pos_in_B = {gene: idx for idx, gene in enumerate(genome2)}
    blocks = []
    i = 0
    while i < len(genome1):
        gene = genome1[i]
        if gene in pos_in_B:
            # potential start of a block
            j = pos_in_B[gene]
            block_len = 1
            # extend block as long if next genes match consecutively
            while i + block_len < len(genome1):
                next_gene = genome1[i + block_len]
                if next_gene in pos_in_B and pos_in_B[next_gene] == j + block_len:
                    block_len += 1
                else:
                    break
            blocks.append(block_len)
            i += block_len  # skip ahead past this block
        else:
            i += 1
    return blocks

#gets number of genomes from the nih database, using its specific format (each data base has to be downloaded, stored in folder; example ATGC0001)
#used in earlier version of the code, not needed anymore; might be useful later
def get_num_genomes():
    colnames = [
        "gene_ID", "genome_ID", "protein_ID", "protein_length",
        "atgc_cog_footprint", "atgc_cog_footprint_length",
        "atgc_cog_ID", "protein_cluster_ID", "match_class"
    ]

    # load the file with correct headers
    df = pd.read_csv("ATGC0001/atgc.cc.csv", names=colnames)

    # drop rows with missing COG annotations
    df = df.dropna(subset=["atgc_cog_ID"])

    # build the presence/absence matrix
    presence_df = df.groupby(["atgc_cog_ID", "genome_ID"]).size().unstack(fill_value=0)
    presence_df = presence_df.applymap(lambda x: 1 if x > 0 else 0)  # Convert counts to 1/0

    # pick a genome column to serve as the root genome
    root_genome_col = presence_df.columns[0]
    root_genome = list(presence_df.index[presence_df[root_genome_col] == 1])
    
    num_genes = len(root_genome)
    return num_genes

#function to load real genome data from file
def get_real_genomes_from_cc():
    #assign column names
    df = pd.read_csv("ATGC0001/atgc.cc.csv", names=[
        "gene_ID", "genome_ID", "protein_ID", "protein_length",
        "atgc_cog_footprint", "atgc_cog_footprint_length",
        "atgc_cog_ID", "protein_cluster_ID", "match_class"
    ])
    #remove rows missing information
    df = df.dropna(subset=["atgc_cog_ID"])
    
    # reconstruct gene order per genome, returns dictionary of genome and its genes
    genome_to_genes = df.groupby("genome_ID")["atgc_cog_ID"].apply(list).to_dict()
    return genome_to_genes

#get root genome
real_genomes = get_real_genomes_from_cc()
root_genome_col = list(real_genomes.keys())[0]
root_genome = real_genomes[root_genome_col]

# testing; using an inversion rate of 0 per gene and 0.1 for inversion/gain per gene for the ATGC0001 database

tree = Phylo.read("ATGC0001/atgc.iq.r.tre", "newick")
sim_genomes = evolve_genome(tree, root_genome, per_gene_gain_rate=0.1, per_gene_loss_rate=0.1, per_gene_inv_rate=0.0)


#comparing two simulated leafs under ATGC0001 structure
leaf1, leaf2 = tree.get_terminals()[0], tree.get_terminals()[1]
sim_blocks = synteny_blocks(sim_genomes[leaf1], sim_genomes[leaf2])
sim_dist = {L: sim_blocks.count(L) for L in set(sim_blocks)}

print("Simulated synteny blocks between", leaf1.name, "and", leaf2.name)
print("Block length distribution:", sim_dist)

#ATGC0001 synteny compariosn of two leafs (as genome lists from real data)
real_genome_ids = list(real_genomes.keys())
real1, real2 = real_genome_ids[0], real_genome_ids[1]
real_blocks = synteny_blocks(real_genomes[real1], real_genomes[real2])
real_dist = {L: real_blocks.count(L) for L in set(real_blocks)}

print("\nReal synteny blocks between", real1, "and", real2)
print("Block length distribution:", real_dist)

# Comparing a real leaf with a simulated one with the data matching the genome_ID
matched_leaf = next((leaf for leaf in tree.get_terminals() if leaf.name in real_genomes), None)

if matched_leaf:
    real = real_genomes[matched_leaf.name]
    simulated = sim_genomes[matched_leaf]
    blocks = synteny_blocks(real, simulated)
    
    print(f"\nComparing Real Genome {matched_leaf.name} with Simulated {matched_leaf.name}:")
    print(f"Synteny block lengths: {blocks}")
    print(f"Block length distribution: { {L: blocks.count(L) for L in set(blocks)} }")
else:
    print("No matching leaf found between real and simulated genome names.")

# fix synteny block analaysis
# use smaller dataset atgc70 for testing its smaller 
# two sep process; read genome data, calc paramaters, analzye syntenyy blocks, real data; the other is simulting the tree, using just tree and root
# distribution of sytnenty blocks of all genome pairs in real data; easier to do for simulted data
# third script: real and simulated genome comparions, it produces wassterstein distance for example
# optimize gain and loss rates, keep inversion rate at 0




