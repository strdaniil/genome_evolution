from synteny_tools import findSyntenyReal2
import re
from typing import Dict, List, Iterable, Tuple

# ---------- parsers for COG/CLS and PTY ----------

# Accepts cog.001562 / COG001562 / 1562  (and tolerates cls.001562 if needed)
_cog_like_re = re.compile(r'(?i)^(?:cog|cls)?\.?0*?(\d+)$')

def _cog_to_int(s: str) -> int:
    s = str(s).strip()
    m = _cog_like_re.match(s)
    if not m:
        raise ValueError(f"Unrecognized COG/CLS token: {s}")
    return int(m.group(1))

def load_genomes_from_pty(
    pty_path: str,
    *,
    prefer_cog: bool = True,
) -> Dict[str, Dict[str, List[int]]]:
    """
    Read an ATGC .pty file and return:
        genomes[genome_ID][partition_ID] -> list[int]   (genes in genomic order)

    Lines look like:
      # 1..2506301 + GCF_000211375.1 NC_015428.1
      NC_015428.1_1|LBUC_RS00015  369..1718  +  GCF_000211375.1  NC_015428.1  WP_...  WP_...  cls.000858  cog.000997
      ...

    We take the last two tokens as cls.* and cog.* (order in file matches your example).
    If prefer_cog=True and cog.* missing/blank, we fall back to cls.*.
    """
    genomes: Dict[str, Dict[str, List[int]]] = {}
    with open(pty_path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                # header lines carry genome/partition span info; we don't need it here
                continue

            parts = line.split()
            # We expect at least:
            # [gene_label, coords, strand, genome_ID, partition_ID, ..., cls_token, cog_token]
            if len(parts) < 8:
                continue  # skip malformed

            genome_id = parts[3]
            part_id   = parts[4]

            # robustly pick CLS/COG tokens from tail
            cls_tok, cog_tok = None, None
            for tok in reversed(parts):
                if tok.lower().startswith("cog."):
                    cog_tok = tok
                    break
            for tok in reversed(parts):
                if tok.lower().startswith("cls."):
                    cls_tok = tok
                    break

            chosen = None
            if prefer_cog and cog_tok:
                chosen = cog_tok
            elif cls_tok:
                chosen = cls_tok
            elif cog_tok:  # if prefer_cog=False or only cog present
                chosen = cog_tok

            if not chosen:
                # no usable label; skip this gene
                continue

            try:
                gene_int = _cog_to_int(chosen)
            except ValueError:
                continue

            genomes.setdefault(genome_id, {}).setdefault(part_id, []).append(gene_int)

    return genomes

# ---------- small accessor & pretty-runner ----------

def get_partition_genes(
    genomes: Dict[str, Dict[str, List[int]]],
    genome_id: str,
    part_id: str
) -> List[int]:
    return list(map(int, genomes[genome_id][part_id]))

def synteny_for_specific_pair(
    genomes: Dict[str, Dict[str, List[int]]],
    g1: str, p1: str,
    g2: str, p2: str
) -> List[List[int]]:
    """
    Prints:
      <g1> <g2> <p1> <p2> <len1> <len2>
      <s1> <e1> <s2> <e2> <L>
      ...
    Returns list [[s1,e1,s2,e2,L], ...] sorted by L desc.
    """
    genome1 = get_partition_genes(genomes, g1, p1)
    genome2 = get_partition_genes(genomes, g2, p2)

    blocks = findSyntenyReal2(genome1[:], genome2[:])
    blocks.sort(key=lambda b: b[4], reverse=True)

    print(f"{g1}\t{g2}\t{p1}\t{p2}\t{len(genome1)}\t{len(genome2)}")
    rows: List[List[int]] = []
    for s1, e1, s2, e2, L in blocks:
        rec = [int(s1), int(e1), int(s2), int(e2), int(L)]
        rows.append(rec)
        print("\t".join(map(str, rec)))
    return rows

# ---------- example call ----------
pty_path = "C:\\Users\\danii\\OneDrive\\Documents\\GitHub\\bio_research\\ATGC0064\\atgc.pty"
genomes = load_genomes_from_pty(pty_path, prefer_cog=True)
# rows = synteny_for_specific_pair(
#     genomes,
#     g1="GCF_000211375.1", p1="NC_015428.1",
#     g2="GCF_000298115.2", p2="NC_018610.1"
# )

# ---------- batch runner over all genomes/partitions ----------

from itertools import combinations, product
from typing import Dict, List, Tuple, Iterable, Optional

def _partition_pairs_all(parts1: Dict[str, List[int]],
                         parts2: Dict[str, List[int]]) -> Iterable[Tuple[str, str]]:
    """
    Generate every (p1, p2) pair. Order by:
      1) max(len1,len2) desc
      2) len1+len2 desc
      3) (p1,p2) tiebreak
    """
    scored = []
    for p1, p2 in product(sorted(parts1.keys()), sorted(parts2.keys())):
        n1, n2 = len(parts1[p1]), len(parts2[p2])
        scored.append((p1, p2, max(n1, n2), n1 + n2))
    scored.sort(key=lambda t: (t[2], t[3], t[0], t[1]), reverse=True)
    for p1, p2, _, __ in scored:
        yield p1, p2

def write_all_pairs_to_file(
    genomes: Dict[str, Dict[str, List[int]]],
    out_path: str,
    *,
    min_block_len: int = 1,
    include_empty: bool = True,   # write header even if no blocks pass the filter
    flush_every: int = 100,       # flush periodically for long runs
) -> None:
    """
    For every genome pair and every partition pair, run findSyntenyReal2
    and write results to `out_path` in your exact format:

        GCF_x  GCF_y  NC_a  NC_b  size_a  size_b
        s1 e1 s2 e2 L
        ...
        <blank line>

    No pairs are skipped; if sequences are empty/unavailable, we still write the header
    if include_empty=True (with zero lengths) followed by a blank line.
    """
    genome_ids = sorted(genomes.keys())
    wrote = 0

    with open(out_path, "w", encoding="utf-8") as fh:
        for g1, g2 in combinations(genome_ids, 2):
            parts1 = genomes.get(g1, {})
            parts2 = genomes.get(g2, {})
            # if a genome somehow has no partitions, still emit nothing for it
            for p1, p2 in _partition_pairs_all(parts1, parts2):
                seq1 = parts1.get(p1, [])
                seq2 = parts2.get(p2, [])

                # header first (always), lengths taken from raw sequences
                fh.write(f"{g1}\t{g2}\t{p1}\t{p2}\t{len(seq1)}\t{len(seq2)}\n")

                rows_written = 0
                if seq1 and seq2:
                    # run synteny on copies (function mutates)
                    blocks = findSyntenyReal2(seq1[:], seq2[:])
                    blocks.sort(key=lambda b: b[4], reverse=True)

                    for s1, e1, s2, e2, L in blocks:
                        L = int(L)
                        if L >= min_block_len:
                            fh.write(f"{int(s1)}\t{int(e1)}\t{int(s2)}\t{int(e2)}\t{L}\n")
                            rows_written += 1

                # blank line separator between comparisons
                fh.write("\n")

                wrote += 1
                if flush_every and (wrote % flush_every == 0):
                    fh.flush()

    print(f"Done. Wrote {wrote} comparison headers to: {out_path}")

# 1) load genomes from your PTY:
pty_path = r"C:\Users\danii\OneDrive\Documents\GitHub\bio_research\ATGC0064\atgc.pty"
genomes = load_genomes_from_pty(pty_path, prefer_cog=True)

# 2) write all genome/partition comparisons:
out_path = r"C:\Users\danii\OneDrive\Documents\GitHub\bio_research\synt.ALL.from_pty.tab"
write_all_pairs_to_file(genomes, out_path, min_block_len=1)

