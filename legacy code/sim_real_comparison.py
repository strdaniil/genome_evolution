from scipy.stats import wasserstein_distance
from simulation import run_simulation, get_real_genomes_from_cc
from calc_real_data import get_pair_blocks
import numpy as np
from Bio import Phylo
from skopt import gp_minimize
from skopt.space import Real
from skopt.utils import use_named_args

DATASET = "ATGC0070"
TREE_FILE = f"{DATASET}/atgc.iq.r.tre"
CC_FILE = f"{DATASET}/atgc.cc.csv"

# Load tree
tree = Phylo.read(TREE_FILE, "newick")



# Load root genome (just use first genome from file)
real_genomes = get_real_genomes_from_cc(CC_FILE)

with open("ATGC0070/atgc.info.tab", "r") as f:
    first_line = f.readline().strip()
    first_value = first_line.split()[0]
root_genome = real_genomes[first_value]

leaves = tree.get_terminals()

best_distance = float('inf')
best_params = None


real_pairs = get_pair_blocks()

sim_pairs = run_simulation(tree, root_genome, 0.1, 0.1, 0.0)

#print(real_pairs, sim_pairs)


def find_error(realpairs, simpairs):
    total_distance = 0
    for pair in realpairs:
        if pair in simpairs:
            real_blocks = realpairs[pair]
            sim_blocks = simpairs[pair]
            dist = wasserstein_distance(real_blocks, sim_blocks)
            total_distance += dist
        else:
            print(f"Missing simulated data for pair: {pair}")
    return total_distance

#print(f"Trying gain={0.1:.3f}, loss={0.1:.3f} -> Distance={find_error(real_pairs, sim_pairs):.4f}")

# optional plot
# all_real_lengths = [length for blocks in real_pairs.values() for length in blocks]
# all_sim_lengths  = [length for blocks in sim_pairs.values() for length in blocks]

# import matplotlib.pyplot as plt
# plt.figure(figsize=(12, 5))

# plt.subplot(1, 2, 1)
# plt.hist(all_real_lengths, bins=range(1, max(all_real_lengths)+2), color='blue', alpha=0.7, edgecolor='black')
# plt.title("Real Synteny Block Lengths")
# plt.xlabel("Block Length")
# plt.ylabel("Frequency")

# plt.subplot(1, 2, 2)
# plt.hist(all_sim_lengths, bins=range(1, max(all_sim_lengths)+2), color='orange', alpha=0.7, edgecolor='black')
# plt.title("Simulated Synteny Block Lengths")
# plt.xlabel("Block Length")
# plt.ylabel("Frequency")

# plt.tight_layout()
# plt.show()


#this is now the training part; will need to change, slow right now

# Define the search space
# search_space = [
#     Real(0.01, 0.5, name="gain"),
#     Real(0.01, 0.5, name="loss")
# ]

# best_distance = float("inf")
# best_params = None

# @use_named_args(search_space)
# def objective(**params):
#     gain = params["gain"]
#     loss = params["loss"]
    
#     sim_pairs = run_simulation(tree, root_genome, gain, loss, inv_rate=0.0)
#     dist = find_error(sim_pairs, real_pairs)
#     print(f"Trying gain={gain:.3f}, loss={loss:.3f} -> Distance={dist:.4f}")
    
#     global best_distance, best_params
#     if dist < best_distance:
#         best_distance = dist
#         best_params = (gain, loss)
    
#     return dist

# #Bayesian optimization
# result = gp_minimize(
#     func=objective,
#     dimensions=search_space,
#     acq_func="EI",  # Expected Improvement
#     n_calls=30,
#     n_initial_points=5,
#     random_state=42
# )

# print("Best parameters:", best_params, "Wasserstein distance:", best_distance)

n = 4 #num of times to repeat sim, right now small num due to computational constraints

# -*- coding: utf-8 -*-
"""
Pooled-Wasserstein pipeline:
- Load tree and real genomes
- Compute real synteny-block length lists per cognate pair
- Run simulation n times; pool simulated block-length counts per pair across runs
- Normalize (pmf) for both real and pooled-sim
- Compute one Wasserstein-1 per pair and sum across pairs

This matches the "pool counts first, then normalize once" instruction.
"""

from collections import Counter, defaultdict
from typing import Dict, List, Tuple
import numpy as np
from scipy.stats import wasserstein_distance
from Bio import Phylo

# --- your modules ---
from simulation import run_simulation, get_real_genomes_from_cc
from calc_real_data import get_pair_blocks

# --------------------
# Config / dataset paths
# --------------------
DATASET = "ATGC0070"
TREE_FILE = f"{DATASET}/atgc.iq.r.tre"
CC_FILE = f"{DATASET}/atgc.cc.csv"
INFO_TAB = f"{DATASET}/atgc.info.tab"

# --------------------
# Helpers
# --------------------
def lengths_to_pmf(lengths: List[int]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a list of block lengths into (values, probs) arrays for a discrete pmf.
    Returns empty arrays if there are no blocks.
    """
    c = Counter(lengths)
    total = sum(c.values())
    if total == 0:
        return np.array([]), np.array([])
    vals = np.fromiter(c.keys(), dtype=int)
    probs = np.fromiter((v / total for v in c.values()), dtype=float)
    return vals, probs

def pooled_counts_to_pmf(counts: Counter) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a Counter(length -> count) into (values, probs) pmf arrays.
    Returns empty arrays if the counter is empty.
    """
    total = sum(counts.values())
    if total == 0:
        return np.array([]), np.array([])
    vals = np.fromiter(counts.keys(), dtype=int)
    probs = np.fromiter((v / total for v in counts.values()), dtype=float)
    return vals, probs

# --------------------
# Load data
# --------------------
tree = Phylo.read(TREE_FILE, "newick")

# Load real genomes
real_genomes = get_real_genomes_from_cc(CC_FILE)

# Choose a root genome (your current approach: first genome id in info.tab)
with open(INFO_TAB, "r") as f:
    first_line = f.readline().strip()
    first_value = first_line.split()[0]
root_genome = real_genomes[first_value]

# All real cognate pairs -> list of real SBL lengths
real_pairs: Dict[Tuple[str, str], List[int]] = get_pair_blocks()

# Precompute real pmfs per pair
real_pmfs: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
for pair, lengths in real_pairs.items():
    real_pmfs[pair] = lengths_to_pmf(lengths)

# --------------------
# Run simulations & pool counts across runs (Mr. Wolf's method)
# --------------------
n_runs = 4  # adjust when you have more compute

# Adjust run_simulation signature as needed; currently: (tree, root_genome, gain, loss, inv_rate)
gain_rate = 0.5
loss_rate = 0.5
inv_rate = 0.0
# If/when you add translocation & block-sized events, extend run_simulation(...) accordingly.

# pair -> Counter(length -> pooled count across runs)
pooled_sim_counts: Dict[Tuple[str, str], Counter] = defaultdict(Counter)

for i in range(n_runs):
    sim_pairs: Dict[Tuple[str, str], List[int]] = run_simulation(
        tree, root_genome, gain_rate, loss_rate, inv_rate
    )
    # Pool raw counts (do NOT normalize per run)
    for pair, sim_lengths in sim_pairs.items():
        pooled_sim_counts[pair].update(sim_lengths)
    print(f"Completed simulation {i+1}/{n_runs}")

# --------------------
# Compute one Wasserstein-1 per pair and sum across pairs
# --------------------

total_wasserstein = 0.0
num_compared_pairs = 0
num_skipped_empty_real = 0
num_skipped_empty_sim = 0

for pair, (real_vals, real_probs) in real_pmfs.items():
    # Skip pairs with no real blocks (define your policy if needed)
    if real_vals.size == 0:
        num_skipped_empty_real += 1
        continue

    sim_counts = pooled_sim_counts.get(pair, None)
    if sim_counts is None or sum(sim_counts.values()) == 0:
        num_skipped_empty_sim += 1
        continue

    sim_vals, sim_probs = pooled_counts_to_pmf(sim_counts)

    # Compute W1 between two discrete pmfs (supports need not match)
    dist = wasserstein_distance(real_vals, sim_vals, u_weights=real_probs, v_weights=sim_probs)
    total_wasserstein += dist
    num_compared_pairs += 1

print(f"Compared pairs: {num_compared_pairs}")
print(f"Skipped (no real blocks): {num_skipped_empty_real}")
print(f"Skipped (no pooled sim blocks): {num_skipped_empty_sim}")
print(f"Sum of Wasserstein distances across pairs (pooled over {n_runs} runs): {total_wasserstein:.6f}")

# --------------------
# (Optional) average distance per compared pair
# --------------------
if num_compared_pairs > 0:
    avg_wasserstein = total_wasserstein / num_compared_pairs
    print(f"Average Wasserstein distance per pair: {avg_wasserstein:.6f}")
else:
    print("No pairs compared; check data generation.")

end_time = time.time()

