#!/usr/bin/env python3
"""
Aggregate evaluate_liftover.py results over simulation replicates (seeds).

usage:
  aggregate_seeds_v2.py BENCH/scen*_seed*/eval --out BENCH/summary_over_seeds_v2.tsv \
                        [--ref-tool snpfix]

Each eval dir must sit in a folder named scen<X>_seed<N> and contain
per_marker.tsv (patched evaluate_liftover.py, with event_id).

Outcomes per marker (add up to 1 in every category, gaps/deletions included):
  correct     correct, or not placed when the marker is unplaceable
  not_placed  no position for a marker that exists in the new assembly
  wrong       wrong position/chromosome, all hits wrong, or any position for an
              unplaceable marker (false placement)
  multi       several positions, the true one among them

Outputs
  --out                 one row per scenario x marker table x category x tool:
                        n markers, n events, n replicates, mean +/- SD over
                        replicates, pooled fraction with 95% Wilson interval
  <out>_paired.tsv      ref tool vs each other tool on the SAME markers:
                        marker level (exact McNemar test) and event level
                        (sign test on events where one tool did better), the
                        event-level test avoids counting markers of one event
                        as independent
Marker tables: "by_category" = background + event_inside markers (Fig. 3),
"realistic" = background markers only (headline accuracy).
"""
import argparse
import math
import os
import re
from collections import Counter, defaultdict

OUTCOME = {"correct": "correct", "correct_unmapped": "correct",
           "unmapped": "not_placed",
           "wrong_pos": "wrong", "wrong_chr": "wrong", "multi_wrong": "wrong",
           "false_placement": "wrong",
           "multi_incl_correct": "multi"}
KINDS = ["correct", "not_placed", "wrong", "multi"]


def wilson(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def binom_two_sided(k, n):
    """Exact two-sided binomial test, p = 0.5."""
    if n == 0:
        return float("nan")
    pk = [math.comb(n, i) * 0.5 ** n for i in range(n + 1)]
    obs = pk[k]
    return min(1.0, sum(p for p in pk if p <= obs * (1 + 1e-9)))


def mean_sd(v):
    if not v:
        return float("nan"), float("nan")
    m = sum(v) / len(v)
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0
    return m, sd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("evaldirs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref-tool", default="snpfix")
    ap.add_argument("--min-markers", type=int, default=5,
                    help="categories with fewer pooled markers are not printed")
    a = ap.parse_args()

    # rows[(scen, table, cat, tool)] -> list of (rep, marker, event, outcome)
    rows = defaultdict(list)
    tools = []
    for d in a.evaldirs:
        m = re.search(r"scen([^_/]+)_seed([^/]+)", os.path.abspath(d))
        if not m:
            raise SystemExit(f"cannot read scenario/seed from path: {d}")
        scen, rep = m.group(1), m.group(2)
        with open(os.path.join(d, "per_marker.tsv")) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if "event_id" not in hdr:
                raise SystemExit(f"{d}/per_marker.tsv has no event_id: rerun the patched "
                                 f"evaluate_liftover.py")
            for line in fh:
                r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                tool = r["tool"]
                if tool not in tools:
                    tools.append(tool)
                out = OUTCOME[r["status"]]
                ev = f"{rep}:{r['event_id']}" if r["event_id"] != "." else f"{rep}:{r['MarkerID']}"
                key = (rep, r["MarkerID"], ev, out)
                if r["marker_set"] in ("background", "event_inside", "NA"):
                    rows[(scen, "by_category", r["category"], tool)].append(key)
                if r["marker_set"] in ("background", "NA"):
                    rows[(scen, "realistic", r["category"], tool)].append(key)
                    rows[(scen, "realistic", "ALL", tool)].append(key)

    # ---------------- summary table
    out = open(a.out, "w")
    out.write("scenario\ttable\tcategory\ttool\tn_markers\tn_events\tn_replicates\t" +
              "\t".join(f"{k}_mean\t{k}_sd" for k in KINDS) +
              "\tcorrect_pooled\tcorrect_ci_low\tcorrect_ci_high\twrong_pooled"
              "\twrong_ci_low\twrong_ci_high\n")
    groups = sorted({k[:3] for k in rows})
    for scen, table, cat in groups:
        n0 = len(rows[(scen, table, cat, tools[0])])
        reps0 = sorted({x[0] for x in rows[(scen, table, cat, tools[0])]})
        ev0 = {x[2] for x in rows[(scen, table, cat, tools[0])] if not x[2].endswith(":.")}
        show = n0 >= a.min_markers
        if show:
            print(f"\nscenario {scen} - {table} - {cat} "
                  f"({n0} markers, {len(ev0)} events/markers as units, {len(reps0)} replicates)")
            print(f"{'tool':12s}" + "".join(f"{k:>18s}" for k in KINDS) +
                  f"{'correct pooled [95% CI]':>30s}")
        for tool in tools:
            lst = rows[(scen, table, cat, tool)]
            per_rep = defaultdict(Counter)
            for rep, _, _, o in lst:
                per_rep[rep][o] += 1
            stats = {}
            for k in KINDS:
                vals = [c[k] / sum(c.values()) for c in per_rep.values() if sum(c.values())]
                stats[k] = mean_sd(vals)
            tot = Counter(o for *_, o in lst)
            n = len(lst)
            cl, ch = wilson(tot["correct"], n)
            wl, wh = wilson(tot["wrong"], n)
            out.write("\t".join(map(str, [scen, table, cat, tool, n,
                      len({x[2] for x in lst}), len(per_rep)])) + "\t" +
                      "\t".join(f"{stats[k][0]:.4f}\t{stats[k][1]:.4f}" for k in KINDS) +
                      f"\t{tot['correct'] / n if n else float('nan'):.4f}\t{cl:.4f}\t{ch:.4f}"
                      f"\t{tot['wrong'] / n if n else float('nan'):.4f}\t{wl:.4f}\t{wh:.4f}\n")
            if show:
                print(f"{tool:12s}" + "".join(f"{stats[k][0]:8.3f} ± {stats[k][1]:5.3f}  "
                                               for k in KINDS) +
                      f"   {tot['correct'] / n:.3f} [{cl:.3f}-{ch:.3f}]")
    out.close()

    # ---------------- paired comparisons: ref tool vs others, same markers
    if a.ref_tool not in tools:
        print(f"\n(--ref-tool {a.ref_tool} not found; paired tests skipped)")
        return
    pfile = re.sub(r"\.tsv$", "", a.out) + "_paired.tsv"
    with open(pfile, "w") as fh:
        fh.write("scenario\ttable\tcategory\tref_tool\tother_tool\tn_markers\t"
                 "only_ref_correct\tonly_other_correct\tmcnemar_p\tn_events\t"
                 "events_ref_better\tevents_other_better\tevents_equal\tsign_test_p\n")
        for scen, table, cat in groups:
            ref = {(r, m): (e, o == "correct") for r, m, e, o in rows[(scen, table, cat, a.ref_tool)]}
            for tool in tools:
                if tool == a.ref_tool:
                    continue
                oth = {(r, m): (e, o == "correct") for r, m, e, o in rows[(scen, table, cat, tool)]}
                common = ref.keys() & oth.keys()
                b = sum(1 for k in common if ref[k][1] and not oth[k][1])
                c = sum(1 for k in common if oth[k][1] and not ref[k][1])
                p_m = binom_two_sided(min(b, c), b + c) if b + c else 1.0
                ev = defaultdict(lambda: [0, 0, 0])
                for k in common:
                    e = ref[k][0]
                    ev[e][0] += ref[k][1]; ev[e][1] += oth[k][1]; ev[e][2] += 1
                better = sum(1 for x in ev.values() if x[0] > x[1])
                worse = sum(1 for x in ev.values() if x[1] > x[0])
                p_s = binom_two_sided(min(better, worse), better + worse) if better + worse else 1.0
                fh.write("\t".join(map(str, [scen, table, cat, a.ref_tool, tool, len(common),
                         b, c, f"{p_m:.4g}", len(ev), better, worse,
                         len(ev) - better - worse, f"{p_s:.4g}"])) + "\n")
    print(f"\nSummary: {a.out}\nPaired tests ({a.ref_tool} vs others): {pfile}")


if __name__ == "__main__":
    main()
