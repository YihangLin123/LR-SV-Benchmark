#!/usr/bin/env python3
"""
eval_MIER.py

Calculate MIER from trio VCF using vcfpy, and report by SVTYPE.
This version recognizes SVTYPE values including TRA and maps TRA -> BND.
"""
import vcfpy
import argparse
from datetime import datetime
from collections import defaultdict
import sys

def now_ts():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[:-3]

def norm_svtype_raw(sv):
    """
    Normalize raw SVTYPE from INFO: handle lists and common variants.
    Map TRA -> BND (treat translocations as breakends).
    """
    if sv is None:
        return None
    if isinstance(sv, (list, tuple)):
        sv = sv[0]
    svs = str(sv).upper()
    if svs == "TRA":
        return "BND"   # map TRA to BND
    return svs

def detect_svtype(record):
    """
    Extract SVTYPE robustly:
      1) try INFO.SVTYPE
      2) fallback to ALT heuristics (breakend symbols or textual clues)
      3) return 'UNK' if cannot identify
    """
    sv = record.INFO.get("SVTYPE")
    svn = norm_svtype_raw(sv)
    if svn:
        return svn

    # fallback: check ALT for breakend symbols or textual hints
    try:
        if record.ALT and len(record.ALT) > 0:
            alt0 = record.ALT[0]
            sval = str(alt0).upper()
            if "[" in sval or "]" in sval:
                return "BND"
            if "DEL" in sval:
                return "DEL"
            if "INS" in sval:
                return "INS"
            if "DUP" in sval:
                return "DUP"
            if "INV" in sval:
                return "INV"
    except Exception:
        pass
    return "UNK"

def norm_gt(gt):
    """ Normalize GT string; convert '|'->'/' and treat '.'/'./.'/'NA' as missing (return None) """
    if gt is None:
        return None
    g = str(gt).strip()
    if g in (".", "./.", ".|.", "", "NA"):
        return None
    return g.replace("|", "/")

def allowed_child(f, m):
    """ Given parental GT strings like '0/1', return set of allowed child GT strings. """
    combos = set()
    afs = f.split("/")
    ams = m.split("/")
    for a in afs:
        for b in ams:
            combos.add("/".join(sorted([a, b])))
    return combos

def main():
    parser = argparse.ArgumentParser(description="Calculate MIER from a trio VCF (vcfpy). Maps TRA -> BND.")
    parser.add_argument("-v", "--vcf", required=True, help="Filtered trio VCF file (last 3 samples are father,mother,child)")
    parser.add_argument("--label", default="F1", help="Sample label to print (default: F1)")
    args = parser.parse_args()

    try:
        reader = vcfpy.Reader.from_path(args.vcf)
    except Exception as e:
        print(f"Failed to open VCF '{args.vcf}': {e}", file=sys.stderr)
        sys.exit(1)

    samples = reader.header.samples.names
    if len(samples) < 3:
        print("VCF has fewer than 3 samples; need a trio.", file=sys.stderr)
        sys.exit(1)

    father_idx = len(samples) - 3
    mother_idx = len(samples) - 2
    child_idx  = len(samples) - 1

    # counters
    total_by_type = defaultdict(int)
    error_by_type = defaultdict(int)
    total_trio_sites = 0
    mendelian_errors = 0

    for rec in reader:
        # get SVTYPE
        svtype = detect_svtype(rec)
        # map uncommon synonyms: treat TRA as BND already in detect_svtype

        # get GTs
        # ensure calls length
        if len(rec.calls) <= child_idx:
            continue
        f_gt = norm_gt(rec.calls[father_idx].data.get("GT"))
        m_gt = norm_gt(rec.calls[mother_idx].data.get("GT"))
        c_gt = norm_gt(rec.calls[child_idx].data.get("GT"))

        # skip missing GTs
        if f_gt is None or m_gt is None or c_gt is None:
            continue

        # count totals
        total_trio_sites += 1
        total_by_type[svtype] += 1

        allowed = allowed_child(f_gt, m_gt)
        if c_gt not in allowed:
            mendelian_errors += 1
            error_by_type[svtype] += 1

    # prepare output ordering
    sv_order = ["DEL","INS","DUP","INV","BND"]
    ts = now_ts()
    label = args.label

    # Print revised table only (you asked to remove the earlier Total/Mendelian/MIER lines)
    print(f"{ts} [INFO] Sample   SV         Total      ME        MIER(%)")
    print(f"{ts} [INFO] ---------------------------------------------")
    for sv in sv_order:
        tot = total_by_type.get(sv, 0)
        err = error_by_type.get(sv, 0)
        mdr = (err / tot * 100.0) if tot>0 else 0.0
        print(f"{ts} [INFO] {label:<8} {sv:<7} {tot:>10} {err:>9} {mdr:12.2f}")
    total_all = sum(total_by_type.values())
    total_err = sum(error_by_type.values())
    total_mdr = (total_err / total_all * 100.0) if total_all>0 else 0.0
    print(f"{ts} [INFO] {label:<8} {'Total':<7} {total_all:>10} {total_err:>9} {total_mdr:12.2f}")
    print(f"{ts} [INFO] Finished.")

if __name__ == "__main__":
    main()
