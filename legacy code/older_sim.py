
def simulate_and_score(
    atgc_dir: str,
    tree_path: str,
    cc_path: str,
    root_mode: str,
    per_gene_gain: float,
    per_gene_loss: float,
    per_gene_inv: float,
    per_gene_trans: float,
    exp_gain: float,
    exp_loss: float,
    exp_inv: float,
    exp_trans: float,
    n_runs: int,
    seed: int | None = None,
    synteny_finder=findSyntenyReal2,
):
    """
    Runs n_runs simulations on the given tree/root genome, pools simulated SBL counts per pair,
    builds pmfs once, computes Wasserstein-1 vs real pmfs for each pair, and sums across pairs.
    """
    # --- Real side: cache to avoid rebuilding for every grid point ---
    cache_key = (tree_path, cc_path, getattr(synteny_finder, "__name__", "finder"))
    try:
        _RPMFS_CACHE  # type: ignore
    except NameError:
        globals()["_RPMFS_CACHE"] = {}

    if cache_key in _RPMFS_CACHE:
        real_pmfs, tree, real_genomes, med_path = _RPMFS_CACHE[cache_key]
    else:
        res = _build_real_pmfs(tree_path, cc_path, synteny_finder=synteny_finder)
        # Support both old (2-return) and new (4-return) versions
        if isinstance(res, tuple) and len(res) == 4:
            real_pmfs, tree, genomes_numeric, med_path = res
            real_genomes = genomes_numeric  # already numeric
        else:
            real_pmfs, tree = res
            real_genomes = _load_real_genomes_from_cc(cc_path)
            med_path = _median_root_to_leaf(tree)
        _RPMFS_CACHE[cache_key] = (real_pmfs, tree, real_genomes, med_path)

    # Root genome (uses cached real_genomes)
    root_genome = _make_root_genome(root_mode, tree, cc_path, real_genomes=real_genomes)

    # Pool simulated counts across runs
    pooled_sim_counts: dict[tuple[str, str], Counter] = defaultdict(Counter)

    # Run simulations
    for _ in range(n_runs):
        sim_pairs = run_simulation(
            tree,
            root_genome,
            per_gene_gain_rate=per_gene_gain,
            per_gene_loss_rate=per_gene_loss,
            per_gene_inv_rate=per_gene_inv,
            per_gene_trans_rate=per_gene_trans,
            gain_exp=exp_gain,
            loss_exp=exp_loss,
            inv_exp=exp_inv,
            trans_exp=exp_trans,
        )
        # Pool counts; standardize pair key ordering without tuple(sorted(...)) overhead
        for (a, b), lens in sim_pairs.items():
            pair = (a, b) if a <= b else (b, a)
            if lens:
                pooled_sim_counts[pair].update(int(x) for x in lens)

    # Compare
    total_w1 = 0.0
    compared = 0
    skipped_real = 0
    skipped_sim = 0

    for pair, (r_vals, r_probs) in real_pmfs.items():
        if r_vals.size == 0:
            skipped_real += 1
            continue
        sim_counts = pooled_sim_counts.get(pair)
        if not sim_counts:
            skipped_sim += 1
            continue
        s_total = sum(sim_counts.values())
        if s_total == 0:
            skipped_sim += 1
            continue

        s_vals = np.fromiter(sorted(sim_counts.keys()), dtype=int)
        s_probs = np.fromiter((sim_counts[v] / s_total for v in s_vals), dtype=float)

        dist = wasserstein_distance(r_vals, s_vals, u_weights=r_probs, v_weights=s_probs)
        total_w1 += float(dist)
        compared += 1

    avg_w1 = (total_w1 / compared) if compared > 0 else float("nan")
    dataset = os.path.basename(os.path.normpath(atgc_dir))

    return {
        "sum_w1": total_w1,
        "avg_w1": avg_w1,
        "n_pairs": compared,
        "skipped_real": skipped_real,
        "skipped_sim": skipped_sim,
        "dataset": dataset,
        "median_root_to_leaf": med_path,
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "params": {
            "root_mode": root_mode,
            "per_gene_gain": per_gene_gain,
            "per_gene_loss": per_gene_loss,
            "per_gene_inv": per_gene_inv,
            "per_gene_trans": per_gene_trans,
            "exp_gain": exp_gain,
            "exp_loss": exp_loss,
            "exp_inv": exp_inv,
            "exp_trans": exp_trans,
            "n_runs": n_runs,
            "seed": seed,
            "tree_path": tree_path,
            "cc_path": cc_path,
        },
    }