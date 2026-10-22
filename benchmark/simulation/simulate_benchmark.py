#!/usr/bin/env python3
"""
simulate_benchmark.py - controlled test assemblies with chip-like markers.

Base assembly (v1): homeologous chromosomes, e.g. the first 10 Mb of wheat
chr1A, chr1B, chr1D (--ref), or a synthetic homeologous triplet (--synthetic).
A new assembly (v2) is derived from v1 with exactly known changes; the true v2
position of every marker is tracked base by base.

Scenario A  "assembly update / homeologs"
    sequence divergence (SNPs + small indels), assembly gaps (N-runs),
    unplaced scaffolds (segments moved to chrUn_scaf*, random orientation,
    the hole closed by a 100-bp N junction). Homeologs 1A/1B/1D are present
    in the base assembly itself.
Scenario B  "synteny disruption"
    stronger divergence, large and small inversions, translocations between
    chromosomes (same and inverted orientation), an intra-chromosomal
    transposition, a reciprocal arm translocation, deletions.

Marker selection follows SNP-array design practice (checked on v1 only):
  * no N in the design flank (+/- --design-flank bp; 35 bp = Axiom 71-mer)
  * GC content of the flank within --gc-min..--gc-max
  * no low-complexity sequence (homopolymer >= 8 bp, dinucleotide repeat >= 10 bp)
  * not repetitive: no 16-mer of the flank occurs more than --max-seed-count
    times in the genome (Axiom "16-mer" repetitiveness check)
  * unique: no other locus (either strand) with >= --max-other-identity
    identity over the whole flank (seed-and-extend with 16-mer seeds,
    tolerating one indel) - the usual "unique BLAST hit" rule
  * optional: only inside gene features of a GFF3 (--genes), as for arrays
    designed from exome capture (90K, 820K, 35K)
  * minimum spacing between markers (--min-spacing)

Marker sets (column marker_set):
  background    chip-like markers at random allowed positions -> realistic
                headline accuracy
  edge_ladder   chip-like markers at fixed distances (0, 5, 10 ... bp) from
                every edge (breakpoints, gap borders, scaffold ends,
                insertion sites, chromosome ends) -> stress test only
  event_inside  chip-like markers inside each event: the midpoint plus
                --inside-markers random positions that pass the design rules
                (redrawn up to --inside-tries times per wanted marker)

Outputs (--outdir): old_ref.fa, new_ref.fa, markers_old.tsv (MarkerID Chr
Position, 1-based), truth.tsv, events.tsv, design_log.tsv (why candidates
were rejected)

Examples
  python simulate_benchmark.py --ref chr1ABD_10Mb.fa --scenario A --outdir simA
  python simulate_benchmark.py --ref chr1ABD_10Mb.fa --scenario B --outdir simB \\
         --genes IWGSC_v1.1_HC.gff3
"""
import argparse
import bisect
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i

SCENARIOS = {
    "A": dict(snp_rate=0.005, indel_rate=0.0005, indel_max=10,
              gaps="200,1000,5000,20000,1000,5000", unplaced="30000,100000,300000",
              junction_gap=100,
              inv_large="", inv_small="", trans_inter="", trans_inter_inv="",
              trans_intra="", deletion="", reciprocal=0),
    "B": dict(snp_rate=0.01, indel_rate=0.001, indel_max=30,
              gaps="", unplaced="", junction_gap=100,
              inv_large="300000,1000000", inv_small="500,2000,10000",
              trans_inter="20000,200000", trans_inter_inv="50000",
              trans_intra="100000", deletion="2000,20000", reciprocal=1),
}
INTERVAL_EVENTS = ("gap", "unplaced_scaffold", "inv_large", "inv_small", "trans_inter",
                   "trans_inter_inv", "trans_intra", "deletion")
# ------------------------------------------------------------------ basics
def revcomp(s: bytes) -> bytes:
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    seqs[name] = b"".join(buf).upper()
                name, buf = line[1:].split()[0].decode(), []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = b"".join(buf).upper()
    return seqs


def write_fasta(path, seqs, width=60):
    with open(path, "w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width].decode() + "\n")


def sizes(s):
    return [int(float(x)) for x in str(s).split(",") if x.strip()] if s else []


def mutate(seq, marker_offsets, rng, snp_rate, indel_rate, indel_max):
    """SNPs + small indels anywhere except on marker bases.
    Returns the new sequence and the new offsets of the markers (same order)."""
    L = len(seq)
    arr = np.frombuffer(seq, dtype=np.uint8).copy()
    mk = np.array(sorted(marker_offsets), dtype=np.int64)
    n = rng.binomial(L, snp_rate) if L else 0
    if n:
        idx = np.unique(rng.integers(0, L, n))
        idx = idx[~np.isin(idx, mk)]
        codes = CODE[arr[idx]]
        ok = codes < 4
        idx, codes = idx[ok], codes[ok]
        arr[idx] = BASES[(codes + rng.integers(1, 4, len(idx))) % 4]
    n = rng.binomial(L, indel_rate) if L > 2 else 0
    pos = np.unique(rng.integers(1, L - 1, n)) if n else []
    chunks, events, cur = [], [], 0
    for q in pos:
        q = int(q)
        if q < cur:
            continue
        k = int(rng.integers(1, indel_max + 1))
        if rng.random() < 0.5:                      # deletion [q, q+k)
            if q + k > L:
                continue
            j = bisect.bisect_left(mk, q)
            if j < len(mk) and mk[j] < q + k:
                continue                            # never delete a marker base
            chunks.append(arr[cur:q].tobytes()); cur = q + k; events.append((q, -k))
        else:                                       # insertion before q
            chunks.append(arr[cur:q].tobytes())
            chunks.append(BASES[rng.integers(0, 4, k)].tobytes())
            cur = q; events.append((q, k))
    chunks.append(arr[cur:].tobytes())
    ev_q = [e[0] for e in events]
    cum = np.cumsum([0] + [e[1] for e in events])
    new = {int(m): int(m) + int(cum[bisect.bisect_right(ev_q, int(m))]) for m in mk}
    return b"".join(chunks), [new[m] for m in marker_offsets]


# ------------------------------------------------------------------ synthetic homeologs
def synthetic_homeologs(length, hdiv, rng, n_te_fam=10, te_per_mb=40):
    """Ancestral chromosome -> three homeologs (chr1A/1B/1D) diverged by hdiv,
    each with its own TE insertions (from a shared TE library)."""
    anc = BASES[rng.choice(4, length, p=[0.27, 0.23, 0.23, 0.27])].tobytes()
    lib = [BASES[rng.integers(0, 4, int(rng.integers(3000, 9000)))].tobytes()
           for _ in range(n_te_fam)]
    out = {}
    for sub in "ABD":
        s, _ = mutate(anc, [], rng, hdiv, hdiv / 10, 10)
        n_te = int(te_per_mb * length / 1e6)
        cuts = np.sort(rng.integers(0, len(s), n_te))
        parts, prev = [], 0
        for c in cuts:
            parts.append(s[prev:c])
            te, _ = mutate(lib[int(rng.integers(len(lib)))], [], rng, 0.02, 0.002, 5)
            parts.append(te if rng.random() < 0.5 else revcomp(te))
            prev = c
        parts.append(s[prev:])
        out[f"chr1{sub}"] = b"".join(parts)
    return out

# ------------------------------------------------------------------ event placement
class Placer:
    def __init__(self, genome, buffer, rng, max_n=0.05):
        self.g, self.buf, self.rng, self.max_n = genome, buffer, rng, max_n
        self.occ = defaultdict(list)
        self.names = list(genome)
        w = np.array([len(genome[c]) for c in self.names], float)
        self.w = w / w.sum()

    def _free(self, c, s, e):
        lo, hi = s - self.buf, e + self.buf
        if lo < 0 or hi > len(self.g[c]):
            return False
        if any(not (hi <= a or lo >= b) for a, b in self.occ[c]):
            return False
        win = self.g[c][s:e] if e > s else self.g[c][max(0, s - 500):s + 500]
        return win.count(b"N") <= self.max_n * max(len(win), 1)

    def place(self, size, chrom=None, exclude=()):
        for _ in range(20000):
            c = chrom or self.names[self.rng.choice(len(self.names), p=self.w)]
            if c in exclude:
                continue
            L = len(self.g[c])
            if L - size - 2 * self.buf <= 0:
                continue
            s = int(self.rng.integers(self.buf, L - size - self.buf))
            if self._free(c, s, s + size):
                self.occ[c].append((s - self.buf, s + size + self.buf))
                return c, s, s + size
        sys.exit(f"ERROR: cannot place an event of {size} bp - genome too small or "
                 f"crowded; reduce event sizes/counts or --buffer")



# ------------------------------------------------------------------ chip-like marker design
class ChipDesigner:
    """Checks candidate SNP positions against array-design rules on v1."""
    K = 16
    LOWCOMP = re.compile(rb"A{8,}|C{8,}|G{8,}|T{8,}|([ACGT]{2})\1{4,}")

    def __init__(self, genome, flank, gc_min, gc_max, max_seed, max_ident, gap_tol=10):
        k = self.K
        self.genome, self.F = genome, flank
        self.gc_min, self.gc_max = gc_min, gc_max
        self.max_seed, self.max_ident, self.gap_tol = max_seed, max_ident, gap_tol
        names = list(genome)
        self.full = (b"N" * k).join(genome[n] for n in names)
        self.off, o = {}, 0
        for n in names:
            self.off[n] = o
            o += len(genome[n]) + k
        self.arr = np.frombuffer(self.full, dtype=np.uint8)
        codes = CODE[self.arr]
        n = len(codes); m = n - k + 1
        bad = codes > 3
        c = np.where(bad, 0, codes).astype(np.uint32)
        fwd = np.zeros(m, dtype=np.uint32)
        for j in range(k):
            fwd = (fwd << np.uint32(2)) | c[j:j + m]
        rcc = (np.uint32(3) - c)[::-1]
        del c
        rev = np.zeros(m, dtype=np.uint32)
        for j in range(k):
            rev = (rev << np.uint32(2)) | rcc[j:j + m]
        del rcc
        rev = rev[::-1]
        self.fwd = fwd
        self.canon = np.minimum(fwd, rev)
        del rev
        cs = np.concatenate(([0], np.cumsum(bad, dtype=np.int64)))
        valid = (cs[k:k + m] - cs[:m]) == 0
        pos = np.nonzero(valid)[0]
        vals = self.canon[pos]
        order = np.argsort(vals, kind="stable")
        self.sv, self.sp = vals[order], pos[order]

    def _identity(self, starts, qa, Lw):
        """Best identity of query qa against a cluster of nearby diagonals,
        allowing one indel (switch between two diagonals)."""
        eqs = []
        for d in starts:
            if d < 0 or d + Lw > len(self.arr):
                continue
            eqs.append(self.arr[d:d + Lw] == qa)
        if not eqs:
            return 0.0
        best = max(int(e.sum()) for e in eqs)
        if len(eqs) > 1:
            pre = [np.concatenate(([0], np.cumsum(e))) for e in eqs]
            for i in range(len(eqs)):
                for j in range(len(eqs)):
                    if i != j:
                        best = max(best, int(np.max(pre[i] + (pre[j][-1] - pre[j]))))
        return best / Lw

    def check(self, c, p):
        """-> (passed, reason, best identity of another locus)"""
        F, k = self.F, self.K
        L = len(self.genome[c])
        if p - F < 0 or p + F >= L:
            return False, "flank_beyond_sequence_end", float("nan")
        ws = self.off[c] + p - F
        Lw = 2 * F + 1
        win = self.full[ws:ws + Lw]
        if b"N" in win:
            return False, "N_in_flank", float("nan")
        gc = (win.count(b"G") + win.count(b"C")) / Lw
        if not (self.gc_min <= gc <= self.gc_max):
            return False, "gc_content", float("nan")
        if self.LOWCOMP.search(win):
            return False, "low_complexity", float("nan")
        js = np.arange(ws, ws + Lw - k + 1)
        q = self.canon[js]
        lo = np.searchsorted(self.sv, q, "left")
        hi = np.searchsorted(self.sv, q, "right")
        if int((hi - lo).max()) > self.max_seed:
            return False, "repetitive", float("nan")
        hits = {"+": set(), "-": set()}
        for i, j in enumerate(js):
            o = j - ws
            for h in self.sp[lo[i]:hi[i]]:
                h = int(h)
                if self.fwd[h] == self.fwd[j]:
                    strand, start = "+", h - o
                else:
                    strand, start = "-", h - (Lw - k - o)
                if start + Lw > ws and start < ws + Lw:   # overlaps the marker itself
                    continue
                hits[strand].add(start)
        best = 0.0
        qf = np.frombuffer(win, dtype=np.uint8)
        qr = np.frombuffer(revcomp(win), dtype=np.uint8)
        for strand, qa in (("+", qf), ("-", qr)):
            starts = sorted(hits[strand])
            cluster = []
            for s in starts + [None]:
                if s is not None and cluster and s - cluster[-1] <= self.gap_tol:
                    cluster.append(s); continue
                if cluster:
                    best = max(best, self._identity(cluster, qa, Lw))
                cluster = [s] if s is not None else []
        if best >= self.max_ident:
            return False, "not_unique", best
        return True, "ok", best


def copy_class(best_identity):
    if best_identity != best_identity:  # nan
        return "NA"
    if best_identity >= 0.90:
        return "close_copy"
    if best_identity >= 0.80:
        return "distant_copy"
    return "no_similar_copy"


def read_gene_intervals(path, genome, feature):
    iv = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[2] != feature or f[0] not in genome:
                continue
            s, e = int(f[3]) - 1, min(int(f[4]), len(genome[f[0]]))
            if e > s:
                iv.append((f[0], s, e))
    return iv


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--ref", help="FASTA with homeologous chromosomes/regions")
    src.add_argument("--synthetic", type=float, help="length of synthetic homeologs (bp)")
    ap.add_argument("--homeolog-div", type=float, default=0.03)
    ap.add_argument("--scenario", choices=["A", "B"], required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=1)
    # array-design rules
    ap.add_argument("--design-flank", type=int, default=35,
                    help="flank checked on each side (35 = Axiom 71-mer, 60 = Infinium)")
    ap.add_argument("--gc-min", type=float, default=0.30)
    ap.add_argument("--gc-max", type=float, default=0.70)
    ap.add_argument("--max-seed-count", type=int, default=50,
                    help="max genome count of any 16-mer in the flank")
    ap.add_argument("--max-other-identity", type=float, default=0.94,
                    help="reject if another locus reaches this identity over the flank")
    ap.add_argument("--genes", help="GFF3; sample background markers only in these features")
    ap.add_argument("--gene-feature", default="gene")
    ap.add_argument("--min-spacing", type=int, default=1000)
    # marker sets
    ap.add_argument("--n-background", type=int, default=3000)
    ap.add_argument("--edge-distances",
                    default="0,5,10,15,20,25,30,35,40,45,50,60,75,100,150,200,300,500,1000")
    ap.add_argument("--no-edge-ladder", action="store_true")
    ap.add_argument("--inside-markers", type=int, default=5,
                    help="markers per event in addition to the midpoint (accepted, not attempted)")
    ap.add_argument("--inside-tries", type=int, default=200,
                    help="max candidate positions per wanted inside marker")
    ap.add_argument("--buffer", type=int, default=20000)
    for p in SCENARIOS["A"]:
        ap.add_argument("--" + p.replace("_", "-"), default=None)
    a = ap.parse_args()

    P = dict(SCENARIOS[a.scenario])
    for p in P:
        v = getattr(a, p)
        if v is not None:
            P[p] = v
    for p in ("snp_rate", "indel_rate"):
        P[p] = float(P[p])
    for p in ("indel_max", "junction_gap", "reciprocal"):
        P[p] = int(float(P[p]))

    rng = np.random.default_rng(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    genome = read_fasta(a.ref) if a.ref else synthetic_homeologs(int(a.synthetic),
                                                                   a.homeolog_div, rng)
    names = list(genome)
    print("v1: " + ", ".join(f"{n} {len(s):,} bp" for n, s in genome.items()), file=sys.stderr)

    # ---------------- events
    placer = Placer(genome, a.buffer, rng)
    events = []

    def add(kind, size_list):
        for sz in sizes(size_list):
            c, s, e = placer.place(sz)
            ev = dict(id=f"{kind}_{len(events)+1}", type=kind, chr=c, start=s, end=e,
                      size=sz, tchr=".", tpos=-1, strand="+")
            if kind in ("inv_large", "inv_small"):
                ev["strand"] = "-"
            elif kind == "unplaced_scaffold":
                ev["strand"] = "-" if rng.random() < 0.5 else "+"
            elif kind in ("trans_inter", "trans_inter_inv"):
                if len(names) < 2:
                    sys.exit("ERROR: translocations need >= 2 chromosomes")
                ev["tchr"], ev["tpos"], _ = placer.place(0, exclude=(c,))
                ev["strand"] = "-" if kind == "trans_inter_inv" else "+"
            elif kind == "trans_intra":
                ev["tchr"], ev["tpos"], _ = placer.place(0, chrom=c)
            events.append(ev)

    add("gap", P["gaps"])
    add("unplaced_scaffold", P["unplaced"])
    add("inv_large", P["inv_large"])
    add("inv_small", P["inv_small"])
    add("trans_inter", P["trans_inter"])
    add("trans_inter_inv", P["trans_inter_inv"])
    add("trans_intra", P["trans_intra"])
    add("deletion", P["deletion"])
    used = set()
    for _ in range(P["reciprocal"]):
        free = [c for c in names if c not in used]
        if len(free) < 2:
            break
        ca, cb = rng.choice(free, 2, replace=False)
        used |= {ca, cb}
        _, ba, _ = placer.place(0, chrom=ca)
        _, bb, _ = placer.place(0, chrom=cb)
        events.append(dict(id=f"reciprocal_{len(events)+1}", type="reciprocal", chr=ca,
                           start=ba, end=ba, size=0, tchr=cb, tpos=bb, strand="+"))

    # ---------------- edges
    edges = defaultdict(dict)
    for c in names:
        edges[c][0] = "chrom_end"
        edges[c][len(genome[c])] = "chrom_end"
    for ev in events:
        if ev["type"] == "reciprocal":
            edges[ev["chr"]][ev["start"]] = "reciprocal"
            edges[ev["tchr"]][ev["tpos"]] = "reciprocal"
            continue
        edges[ev["chr"]][ev["start"]] = ev["type"]
        edges[ev["chr"]][ev["end"]] = ev["type"]
        if ev["tpos"] >= 0:
            edges[ev["tchr"]][ev["tpos"]] = "insertion_site"
    cuts = {c: sorted(edges[c]) for c in names}

    # ---------------- chip-like markers
    print("indexing v1 for marker design ...", file=sys.stderr)
    designer = ChipDesigner(genome, a.design_flank, a.gc_min, a.gc_max,
                            a.max_seed_count, a.max_other_identity)
    mpos = defaultdict(set)
    source, best_id = {}, {}
    log = Counter()

    def try_add(c, p, set_name):
        if (c, p) in source:
            return True
        passed, reason, b = designer.check(c, p)
        log[(set_name, reason)] += 1
        if passed:
            mpos[c].add(p); source[(c, p)] = set_name; best_id[(c, p)] = b
        return passed

    # background first (they define the realistic set)
    if a.genes:
        iv = read_gene_intervals(a.genes, genome, a.gene_feature)
        if not iv:
            sys.exit(f"ERROR: no '{a.gene_feature}' features on {names} in {a.genes}")
        w = np.array([e - s for _, s, e in iv], float); w /= w.sum()
        def draw():
            c, s, e = iv[rng.choice(len(iv), p=w)]
            return c, int(rng.integers(s, e))
    else:
        tot = sum(len(genome[c]) for c in names)
        w = np.array([len(genome[c]) / tot for c in names])
        def draw():
            c = names[rng.choice(len(names), p=w)]
            return c, int(rng.integers(0, len(genome[c])))
    accepted = defaultdict(list)
    n_bg = 0
    for _ in range(a.n_background * 200):
        if n_bg >= a.n_background:
            break
        c, p = draw()
        lst = accepted[c]
        j = bisect.bisect_left(lst, p)
        if (j < len(lst) and lst[j] - p < a.min_spacing) or \
           (j > 0 and p - lst[j - 1] < a.min_spacing):
            log[("background", "min_spacing")] += 1
            continue
        if (c, p) not in source and try_add(c, p, "background"):
            bisect.insort(lst, p); n_bg += 1
    if n_bg < a.n_background:
        print(f"WARNING: only {n_bg}/{a.n_background} background markers passed the "
              f"design rules", file=sys.stderr)

    if not a.no_edge_ladder:
        dists = sizes(a.edge_distances)
        for c in names:
            for b in cuts[c]:
                for d in dists:
                    for p in (b - 1 - d, b + d):
                        if 0 <= p < len(genome[c]):
                            try_add(c, p, "edge_ladder")
    # Inside markers: keep drawing random positions until the event has
    # 1 + --inside-markers markers that pass the design rules (the midpoint is
    # tried first). In repeat-rich sequence most candidates fail, so a fixed
    # number of draws would leave many events with few or no markers.
    short_events = 0
    for ev in events:
        if ev["size"] <= 0:
            continue
        want = 1 + a.inside_markers
        got = set()
        cands = [ev["start"] + ev["size"] // 2]
        tries = 0
        while len(got) < want and tries < a.inside_tries * want:
            p = cands.pop() if cands else int(rng.integers(ev["start"], ev["end"]))
            tries += 1
            if p in got:
                continue
            if try_add(ev["chr"], p, "event_inside"):
                got.add(p)
        ev["n_inside"] = len(got)
        if len(got) < want:
            short_events += 1
    if short_events:
        print(f"NOTE: {short_events} events got fewer than {1 + a.inside_markers} inside "
              f"markers (small or repetitive events; see n_inside in events.tsv)",
              file=sys.stderr)

    # ---------------- layout of v2
    layout = {c: [dict(src=(c, s, e), strand="+", kind="seq")
                  for s, e in zip(cuts[c][:-1], cuts[c][1:])] for c in names}

    def find(c, start):
        for nc, lst in layout.items():
            for i, it in enumerate(lst):
                if it["kind"] != "N" and it["src"][0] == c and it["src"][1] == start:
                    return nc, i
        raise KeyError((c, start))

    n_scaf = 0
    for ev in events:
        t = ev["type"]
        if t in ("gap", "deletion", "inv_large", "inv_small", "unplaced_scaffold"):
            nc, i = find(ev["chr"], ev["start"])
            if t == "gap":
                layout[nc][i]["kind"] = "gap"
            elif t == "deletion":
                layout[nc].pop(i)
            elif t in ("inv_large", "inv_small"):
                layout[nc][i]["strand"] = "-"
            else:
                item = layout[nc][i]
                layout[nc][i] = dict(kind="N", length=P["junction_gap"])
                item["strand"] = ev["strand"]
                n_scaf += 1
                ev["tchr"] = f"chrUn_scaf{n_scaf}"
                layout[ev["tchr"]] = [item]
    for ev in events:
        if ev["type"] in ("trans_inter", "trans_inter_inv", "trans_intra"):
            nc, i = find(ev["chr"], ev["start"])
            item = layout[nc].pop(i)
            item["strand"] = ev["strand"]
            tc, j = find(ev["tchr"], ev["tpos"])
            layout[tc].insert(j, item)
    for ev in events:
        if ev["type"] == "reciprocal":
            ca, ia = find(ev["chr"], ev["start"])
            cb, ib = find(ev["tchr"], ev["tpos"])
            ta, tb = layout[ca][ia:], layout[cb][ib:]
            layout[ca] = layout[ca][:ia] + tb
            layout[cb] = layout[cb][:ib] + ta

    # ---------------- build v2 and the truth
    msorted = {c: sorted(v) for c, v in mpos.items()}
    truth = defaultdict(list)
    v2 = {}
    for nc, items in layout.items():
        parts, off = [], 0
        for it in items:
            if it["kind"] == "N":
                parts.append(b"N" * it["length"]); off += it["length"]; continue
            c, s, e = it["src"]
            lo = bisect.bisect_left(msorted.get(c, []), s)
            hi = bisect.bisect_left(msorted.get(c, []), e)
            ms = msorted.get(c, [])[lo:hi]
            if it["kind"] == "gap":
                parts.append(b"N" * (e - s)); off += e - s; continue
            seq, new_off = mutate(genome[c][s:e], [m - s for m in ms], rng,
                                  P["snp_rate"], P["indel_rate"], P["indel_max"])
            if it["strand"] == "-":
                seq = revcomp(seq)
                new_off = [len(seq) - 1 - o for o in new_off]
            for m, o in zip(ms, new_off):
                truth[(c, m)].append((nc, off + o + 1,
                                      "REVERSE" if it["strand"] == "-" else "FORWARD"))
            parts.append(seq); off += len(seq)
        v2[nc] = b"".join(parts)

    # ---------------- annotate and write
    def category(c, p):
        for ev in events:
            if ev["type"] in INTERVAL_EVENTS and ev["chr"] == c and ev["start"] <= p < ev["end"]:
                return ev["type"], ev["id"]
        return "collinear", "."

    def nearest_edge(c, p):
        b = cuts[c]
        j = bisect.bisect_right(b, p)
        cand = []
        if j > 0:
            cand.append((p - b[j - 1], b[j - 1]))
        if j < len(b):
            cand.append((b[j] - 1 - p, b[j]))
        d, cut = min(cand)
        return d, edges[c][cut]

    write_fasta(os.path.join(a.outdir, "old_ref.fa"), genome)
    write_fasta(os.path.join(a.outdir, "new_ref.fa"), v2)
    summary = Counter()
    with open(os.path.join(a.outdir, "markers_old.tsv"), "w") as chip, \
         open(os.path.join(a.outdir, "truth.tsv"), "w") as tr:
        chip.write("MarkerID\tChr\tPosition\n")
        tr.write("MarkerID\told_chr\told_pos\tnew_chr\tnew_pos\tnew_strand\tcategory\t"
                 "event_id\tdist_to_bp\tedge_type\tcopy_class\tbest_other_identity\t"
                 "expected\tmarker_set\n")
        n = 0
        for c in names:
            for p in msorted.get(c, []):
                n += 1
                mid = f"M{n:06d}"
                cat, eid = category(c, p)
                d, etype = nearest_edge(c, p)
                b = best_id[(c, p)]
                t = truth.get((c, p), [])
                expected = "place" if t else "unplaceable"
                summary[(source[(c, p)], cat, expected)] += 1
                chip.write(f"{mid}\t{c}\t{p+1}\n")
                tr.write("\t".join(map(str, [
                    mid, c, p + 1,
                    ";".join(x[0] for x in t) or "NA",
                    ";".join(str(x[1]) for x in t) or "NA",
                    ";".join(x[2] for x in t) or "NA",
                    cat, eid, d, etype, copy_class(b), f"{b:.3f}", expected,
                    source[(c, p)]])) + "\n")
    with open(os.path.join(a.outdir, "events.tsv"), "w") as fh:
        fh.write("event_id\ttype\tchr\tstart0\tend0\tsize\ttarget\ttarget_pos0\tstrand\tn_inside\n")
        for ev in events:
            fh.write("\t".join(map(str, [ev[k] for k in ("id", "type", "chr", "start", "end",
                     "size", "tchr", "tpos", "strand")] + [ev.get("n_inside", 0)])) + "\n")
    with open(os.path.join(a.outdir, "design_log.tsv"), "w") as fh:
        fh.write("marker_set\tresult\tcandidates\n")
        for (s, r), k in sorted(log.items()):
            fh.write(f"{s}\t{r}\t{k}\n")

    print("v2: " + ", ".join(f"{k} {len(s):,} bp" for k, s in v2.items()), file=sys.stderr)
    print(f"{n} markers, {len(events)} events -> {a.outdir}/", file=sys.stderr)
    print("design filter (background candidates): " + ", ".join(
        f"{r} {k}" for (s, r), k in sorted(log.items()) if s == "background"), file=sys.stderr)
    for (ms, cat, ex), k in sorted(summary.items()):
        print(f"  {ms:13s} {cat:18s} {ex:12s} {k}", file=sys.stderr)


if __name__ == "__main__":
    main()
