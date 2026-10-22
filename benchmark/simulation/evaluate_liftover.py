#!/usr/bin/env python3
"""
Score marker-position transfer tools against truth.tsv from simulate_sv.py.

usage:
  evaluate_liftover.py --truth sim/truth.tsv \
      --pred main=out_main/chip_new.txt \
      --pred crossmap=out_cm/chip_crossmap.tsv \
      --pred rustybam=out_rb/chip_rustybam.tsv \
      --pred pyliftover=out_py/chip_pylift.tsv \
      --outdir eval

Prediction files: tab/space separated. If a header is present, columns named
MarkerID/Name/ID, Chr/Chrom/Chromosome and Position/Pos are used; otherwise
columns 1-3. Non-integer positions (n.a., NA, Failed_BLAST ...) = unmapped.
Several rows for one marker = multiple placements.

Per marker and tool the status is one of:
  correct            one placement, right chr and position (within --tol);
                     for duplicated markers: exactly the set of true copies
  multi_incl_correct several placements, one of them right
  multi_wrong        several placements, none right
  wrong_pos          one placement, right chr, wrong position
  wrong_chr          one placement, wrong chromosome
  unmapped           no placement although the marker exists in the new assembly
  correct_unmapped   no placement, marker was deleted (correct)
  false_placement    placement given for a deleted marker

Every summary table also has outcome fractions that add up to 1 for every
category, including gaps and deletions, where the correct answer is
"not placed":
  frac_correct     correct + correct_unmapped
  frac_not_placed  unmapped (marker exists in the new assembly, no position)
  frac_wrong       wrong_pos + wrong_chr + multi_wrong + false_placement
  frac_multi       multi_incl_correct (ambiguous, true locus among the hits)
plus n_placeable / n_unplaceable, so tables are never NaN for gap/deletion.
per_marker.tsv has an event_id column for event-level statistics.
"""
import argparse
import os
from collections import Counter, defaultdict

BINS = [(0, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 74), (75, 99),
        (100, 199), (200, 499), (500, 999), (1000, 10**12)]
ID_COLS = {"markerid", "name", "marker", "id", "snp", "snp_id"}
CHR_COLS = {"chr", "chrom", "chromosome", "seqid"}
POS_COLS = {"position", "pos", "bp"}


def to_int(x):
    try:
        return int(float(x))
    except ValueError:
        return None


def read_pred(path):
    preds = defaultdict(set)
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        return preds
    split = (lambda l: l.split("\t")) if "\t" in lines[0] else (lambda l: l.split())
    first = [x.strip().lower() for x in split(lines[0])]
    ci, cc, cp = 0, 1, 2
    if len(first) > 2 and to_int(first[2]) is None and first[2] not in ("n.a.", "na", "nan"):
        for i, h in enumerate(first):
            if h in ID_COLS: ci = i
            elif h in CHR_COLS: cc = i
            elif h in POS_COLS: cp = i
        lines = lines[1:]
    for l in lines:
        f = split(l)
        if len(f) <= max(ci, cc, cp):
            preds[f[0]]  # touch: marker seen, unmapped
            continue
        pos = to_int(f[cp])
        chrom = f[cc].replace("lcl|", "")
        if pos is None or chrom.lower() in ("n.a.", "na", "*", "."):
            preds[f[ci]]
            continue
        preds[f[ci]].add((chrom, pos))
    return preds


def read_truth(path):
    truth = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for l in fh:
            r = dict(zip(hdr, l.rstrip("\n").split("\t")))
            if r["new_chr"] == "NA":
                r["true"] = []
            else:
                r["true"] = list(zip(r["new_chr"].split(";"),
                                     map(int, r["new_pos"].split(";"))))
            r["dist_to_bp"] = int(r["dist_to_bp"])
            truth[r["MarkerID"]] = r
    return truth


def classify(true, pred, tol):
    if not true:
        return ("correct_unmapped" if not pred else "false_placement"), None
    if not pred:
        return "unmapped", None
    hit = any(pc == tc and abs(pp - tp) <= tol for pc, pp in pred for tc, tp in true)
    if len(pred) > 1:
        m = lambda a, b: a[0] == b[0] and abs(a[1] - b[1]) <= tol
        if len(true) > 1 and all(any(m(p, t) for t in true) for p in pred) \
                and all(any(m(p, t) for p in pred) for t in true):
            return "correct", None  # duplicated marker: all true copies, nothing else
        return ("multi_incl_correct" if hit else "multi_wrong"), None
    if hit:
        return "correct", None
    (pc, pp), = pred
    same = [tp for tc, tp in true if tc == pc]
    if same:
        return "wrong_pos", min((pp - tp for tp in same), key=abs)
    return "wrong_chr", None


def ACC(c):
    """correct answers: exact single placement, or no placement for an
    unplaceable marker (inside a gap or a deletion)."""
    return c["correct"] + c["correct_unmapped"]


WRONG = ("wrong_pos", "wrong_chr", "multi_wrong", "false_placement")
OUTCOME_HDR = "\tn_placeable\tn_unplaceable\tfrac_correct\tfrac_not_placed\tfrac_wrong\tfrac_multi"


def outcome_cols(c):
    """Outcome fractions over ALL markers of a group (placeable or not)."""
    n = sum(c.values())
    unpl = c["correct_unmapped"] + c["false_placement"]
    if not n:
        return "\t0\t0\tnan\tnan\tnan\tnan"
    f = lambda k: f"{k / n:.4f}"
    return (f"\t{n - unpl}\t{unpl}\t{f(ACC(c))}\t{f(c['unmapped'])}"
            f"\t{f(sum(c[x] for x in WRONG))}\t{f(c['multi_incl_correct'])}")


def bin_of(d):
    for lo, hi in BINS:
        if lo <= d <= hi:
            return f"{lo}-{hi}" if hi < 10**12 else f">={lo}"
    return "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--pred", action="append", required=True, help="name=path")
    ap.add_argument("--tol", type=int, default=0, help="position tolerance (bp)")
    ap.add_argument("--outdir", default="eval")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    truth = read_truth(a.truth)
    tools = [p.split("=", 1) for p in a.pred]
    per_marker = []
    offsets = defaultdict(Counter)
    for name, path in tools:
        pred = read_pred(path)
        extra = set(pred) - set(truth)
        if extra:
            print(f"[{name}] {len(extra)} marker IDs not in truth (ignored), e.g. {sorted(extra)[:3]}")
        for mid, r in truth.items():
            st, off = classify(r["true"], pred.get(mid, set()), a.tol)
            if off is not None:
                offsets[(name, r["category"])][off] += 1
            per_marker.append((name, mid, r["category"], r["dist_to_bp"], st,
                               r.get("edge_type", "NA"), r.get("copy_class", "NA"),
                               r.get("marker_set", "NA"), r.get("event_id", ".")))

    with open(os.path.join(a.outdir, "per_marker.tsv"), "w") as fh:
        fh.write("tool\tMarkerID\tcategory\tdist_to_bp\tstatus\tedge_type\tcopy_class\t"
                 "marker_set\tevent_id\n")
        for row in per_marker:
            fh.write("\t".join(map(str, row)) + "\n")

    statuses = ["correct", "multi_incl_correct", "multi_wrong", "wrong_pos", "wrong_chr",
                "unmapped", "correct_unmapped", "false_placement"]
    names = [t[0] for t in tools]

    # ---- by category
    # marker sets: background = realistic chip; event_inside = extra markers in
    # events; edge_ladder = stress test (only used in the distance tables)
    has_sets = any(r[7] != "NA" for r in per_marker)
    realistic = lambda r: (not has_sets) or r[7] == "background"
    by_cat = defaultdict(Counter)
    for r in per_marker:
        if (not has_sets) or r[7] in ("background", "event_inside"):
            by_cat[(r[0], r[2])][r[4]] += 1
    cats = sorted({c for _, c in by_cat})
    with open(os.path.join(a.outdir, "summary_by_category.tsv"), "w") as fh:
        fh.write("tool\tcategory\tn\t" + "\t".join(statuses) + "\taccuracy" + OUTCOME_HDR + "\n")
        print("\nAccuracy by event type - background + event_inside markers "
              "(behaviour per event, NOT overall accuracy)")
        print(f"{'category':18s}" + "".join(f"{n:>12s}" for n in names))
        for cat in cats:
            line = f"{cat:18s}"
            for tool in names:
                c = by_cat[(tool, cat)]
                n = sum(c.values())
                acc = ACC(c) / n if n else float("nan")
                fh.write("\t".join(map(str, [tool, cat, n] + [c[s] for s in statuses]))
                         + f"\t{acc:.4f}" + outcome_cols(c) + "\n")
                line += f"{acc:12.3f}"
            print(line + f"   (n={sum(by_cat[(names[0], cat)].values())})")

    # ---- by distance to breakpoint (markers that exist in the new assembly)
    by_d = defaultdict(Counter)
    for tool, mid, cat, d, st, *_ in per_marker:
        if truth[mid]["true"]:
            by_d[(tool, bin_of(d))][st] += 1
    bins = [bin_of(lo) for lo, _ in BINS]
    with open(os.path.join(a.outdir, "summary_by_distance.tsv"), "w") as fh:
        fh.write("tool\tdist_bin\tn\t" + "\t".join(statuses) + "\taccuracy" + OUTCOME_HDR + "\n")
        # all marker sets: each distance bin is reported on its own, so the
        # stress-test markers do not distort anything here
        print("\nAccuracy vs distance to nearest breakpoint (markers present in new assembly)")
        print(f"{'distance':18s}" + "".join(f"{n:>12s}" for n in names))
        for b in bins:
            line = f"{b:18s}"
            for tool in names:
                c = by_d[(tool, b)]
                n = sum(c.values())
                acc = ACC(c) / n if n else float("nan")
                fh.write("\t".join(map(str, [tool, b, n] + [c[s] for s in statuses]))
                         + f"\t{acc:.4f}" + outcome_cols(c) + "\n")
                line += f"{acc:12.3f}"
            print(line + f"   (n={sum(by_d[(names[0], b)].values())})")

    # ---- by copy class and by edge type (markers within 100 bp of an edge)
    def group_summary(label, key_fn, fname):
        grp = defaultdict(Counter)
        for row in per_marker:
            ks = key_fn(row)
            if ks is None:
                continue
            for k in ([ks] if isinstance(ks, str) else ks):
                grp[(row[0], k)][row[4]] += 1
        keys = sorted({k for _, k in grp})
        if not keys or keys == ["NA"]:
            return
        with open(os.path.join(a.outdir, fname), "w") as fh:
            fh.write("tool\tgroup\tn\t" + "\t".join(statuses) + "\taccuracy" + OUTCOME_HDR + "\n")
            print(f"\nAccuracy by {label}")
            print(f"{'group':18s}" + "".join(f"{n:>12s}" for n in names))
            for k in keys:
                line = f"{k:18s}"
                for tool in names:
                    c = grp[(tool, k)]; n = sum(c.values())
                    acc = ACC(c) / n if n else float("nan")
                    fh.write("\t".join(map(str, [tool, k, n] + [c[x] for x in statuses]))
                             + f"\t{acc:.4f}" + outcome_cols(c) + "\n")
                    line += f"{acc:12.3f}"
                print(line + f"   (n={sum(grp[(names[0], k)].values())})")

    group_summary("category - REALISTIC SET (background markers only, headline numbers)",
                  lambda r: [r[2], "ALL"] if r[7] == "background" else None,
                  "summary_realistic.tsv")
    group_summary("closest other copy in v1 (homeolog similarity) - background markers",
                  lambda r: r[6] if realistic(r) else None, "summary_by_copy_class.tsv")
    group_summary("edge type, markers <=100 bp from the edge (stress test)",
                  lambda r: r[5] if r[3] <= 100 else None, "summary_by_edge_type.tsv")

    # ---- error profile
    real_st = defaultdict(Counter)
    for r in per_marker:
        if realistic(r):
            real_st[r[0]][r[4]] += 1
    print("\nStatus counts - background markers (realistic set)")
    print(f"{'status':20s}" + "".join(f"{n:>12s}" for n in names))
    for s in statuses:
        print(f"{s:20s}" + "".join(f"{real_st[t][s]:12d}" for t in names))
    print("Summary (background markers, markers that exist in the new assembly):")
    for t in names:
        c = real_st[t]
        placeable = sum(c[x] for x in statuses if x not in ("correct_unmapped", "false_placement"))
        wrong = c["wrong_pos"] + c["wrong_chr"] + c["multi_wrong"]
        if placeable:
            print(f"  {t:12s} correct {c['correct']/placeable:.3f}   "
                  f"not placed {c['unmapped']/placeable:.3f}   "
                  f"placed WRONG {wrong/placeable:.3f}   "
                  f"multi incl. true {c['multi_incl_correct']/placeable:.3f}")
        neg = c["correct_unmapped"] + c["false_placement"]
        if neg:
            print(f"  {'':12s} deleted/gap markers wrongly placed: {c['false_placement']}/{neg}")

    print("\nMost common signed offsets for wrong_pos (a constant offset = coordinate bug,")
    print("e.g. +/-1 = 0/1-based mix-up, a value tied to flank size = strand handling):")
    for (tool, cat), cnt in sorted(offsets.items()):
        print(f"  {tool:12s} {cat:16s} " +
              ", ".join(f"{o:+d} (x{k})" for o, k in cnt.most_common(3)))

    # ---- plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(14, 5))
        has_bg = any(r[7] == "background" for r in per_marker)
        bar = defaultdict(Counter)
        for r in per_marker:
            if not has_bg or r[7] == "background":
                bar[(r[0], r[2])][r[4]] += 1
        plot_cats = sorted({k for _, k in bar})
        w = 0.8 / len(names)
        for i, tool in enumerate(names):
            vals = []
            for cat in plot_cats:
                c = bar[(tool, cat)]; n = sum(c.values())
                vals.append(ACC(c) / n if n else 0)
            ax[0].bar([x + i * w for x in range(len(plot_cats))], vals, w, label=tool)
        ax[0].set_xticks([x + w * (len(names) - 1) / 2 for x in range(len(plot_cats))])
        ax[0].set_xticklabels(plot_cats, rotation=35, ha="right")
        ax[0].set_ylabel("fraction correct"); ax[0].set_ylim(0, 1.05)
        ax[0].set_title("Accuracy by SV class" + (" (background markers)" if has_bg else ""))
        ax[0].legend()
        for tool in names:
            ys = []
            for b in bins:
                c = by_d[(tool, b)]; n = sum(c.values())
                ys.append(ACC(c) / n if n else float("nan"))
            ax[1].plot(range(len(bins)), ys, marker="o", label=tool)
        ax[1].set_xticks(range(len(bins))); ax[1].set_xticklabels(bins, rotation=35, ha="right")
        ax[1].set_xlabel("distance to nearest breakpoint (bp)")
        ax[1].set_ylabel("fraction correct"); ax[1].set_ylim(0, 1.05)
        ax[1].set_title("Accuracy vs distance to breakpoint"); ax[1].legend()
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "benchmark.png"), dpi=150)
        print(f"\nPlot: {a.outdir}/benchmark.png")
    except ImportError:
        print("\n(matplotlib not available: no plot)")


if __name__ == "__main__":
    main()
