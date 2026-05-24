#!/usr/bin/env python3
import sys
import argparse
import logging
import time

USAGE = """\
    Evaluate SV callset (Mendelian inconsistency rate) using position/length/breakpoint matching.
"""

def pase_info(seq):
    """
    Parse INFO field (simple parser).
    Returns dict with keys: SVLEN, END, SVTYPE, RE, CHR2 (defaults present).
    """
    info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "RE": 0, "CHR2": ''}
    for kv in seq.split(';'):
        if '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        if k in ["SVLEN", "END", "RE"]:
            try:
                info[k] = abs(int(float(v)))
            except Exception:
                # keep default on parse error
                pass
        elif k == "CHR2":
            info[k] = v
        elif k == "SVTYPE":
            # keep first 3 chars for compatibility
            info[k] = v[0:3]
    return info


def load_callset(path):
    """
    Load a VCF-like callset into a dict:
    - For DEL/INS/DUP/INV: append [chr, pos, END, SVLEN, matched_flag]
    - For BND: append [chr, pos, CHR2, END, form, matched_flag]
    """
    callset = dict()
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line[0] == '#':
                continue
            seq = line.split('\t')
            if len(seq) < 8:
                logging.warning("Skipping malformed line (too few columns): %s", line[:100])
                continue

            chr = seq[0]
            try:
                pos = int(seq[1])
            except Exception:
                logging.warning("Skipping line with non-integer POS: %s", line[:100])
                continue

            info = pase_info(seq[7])

            if info['SVTYPE'] in ['DEL', 'INS', 'DUP', 'INV']:
                if info['SVTYPE'] not in callset:
                    callset[info['SVTYPE']] = list()
                if info['SVLEN'] == 0:
                    # if END is present, approximate SVLEN
                    try:
                        info['SVLEN'] = info['END'] - pos + 1
                        if info['SVLEN'] < 0:
                            info['SVLEN'] = 0
                    except Exception:
                        info['SVLEN'] = 0
                # store: [chr, pos, END, SVLEN, matched_flag]
                callset[info['SVTYPE']].append([chr, pos, info['END'], info['SVLEN'], 0])

            elif info['SVTYPE'] == "BND":
                # parse ALT field (seq[4]) to get chr2 and pos2 and form
                alt = seq[4]
                chr2 = ''
                pos2 = 0
                form = ''
                try:
                    # attempt to parse four common BND forms used previously
                    if alt.startswith(']'):
                        # form like ]chr2:pos]N or ]chr2:pos]...
                        form = ']]N'
                        chr2 = alt.split(':')[0][1:]
                        pos2 = int(alt.split(':')[1].split(']')[0])
                    elif alt.startswith('['):
                        form = '[[N'
                        chr2 = alt.split(':')[0][1:]
                        pos2 = int(alt.split(':')[1].split('[')[0])
                    else:
                        # alt starts with base then bracket e.g. N]chr2:pos] or N[chr2:pos[
                        if len(alt) > 1 and alt[1] == ']':
                            form = 'N]]'
                            chr2 = alt.split(':')[0][2:]
                            pos2 = int(alt.split(':')[1].split(']')[0])
                        else:
                            form = 'N[['
                            chr2 = alt.split(':')[0][2:]
                            pos2 = int(alt.split(':')[1].split('[')[0])
                except Exception:
                    # fallback: try to find "chr:pos" inside alt using simple search
                    try:
                        # extract digits after ':' up to non-digit
                        parts = alt.replace(']', ':').replace('[', ':').split(':')
                        # find last pair that looks like chr,pos
                        if len(parts) >= 2:
                            chr2 = parts[-2]
                            pos2 = int(parts[-1])
                    except Exception:
                        chr2 = ''
                        pos2 = 0

                if info['SVTYPE'] not in callset:
                    callset[info['SVTYPE']] = list()
                if info['END'] == 0:
                    info['CHR2'] = chr2
                    info['END'] = pos2

                # store: [chr, pos, CHR2, END, form, matched_flag]
                callset[info['SVTYPE']].append([chr, pos, info['CHR2'], info['END'], form, 0])

            else:
                # if SVTYPE empty or unrecognized, skip
                continue

    return callset


def eva_record(call_A, call_B, bias, offect):
    """
    Compare call_A vs call_B and mark matches in call_B by setting matched_flag to 1.
    Matching rules:
    - INS: position distance <= offect and length ratio >= bias
    - BND: chr2 equal, form equal, both breakpoints within offect
    - DEL/DUP/INV: intervals overlap (allowing offect extension) and length ratio >= bias
    Note: This function modifies call_B in-place (marking matches).
    """
    # Iterate svtypes present in call_A
    for svtype in call_A:
        if svtype not in call_B:
            continue
        for i in call_B[svtype]:
            # i is an element from call_B[svtype]; structure differs for BND vs others
            for j in call_A[svtype]:
                # require same chromosome for non-BND types
                if svtype == 'BND':
                    # BND layout in call_B: [chr, pos, CHR2, END, form, matched_flag]
                    # Similarly for call_A entries
                    try:
                        if i[2] == j[2] and i[4] == j[4]:
                            if abs(i[1] - j[1]) <= offect and abs(i[3] - j[3]) <= offect:
                                i[-1] = 1
                    except Exception:
                        # be permissive on malformed entries
                        continue
                elif svtype == 'INS':
                    # layout: [chr, pos, END, SVLEN, matched_flag]
                    try:
                        if i[0] != j[0]:
                            continue
                        if abs(i[1] - j[1]) <= offect:
                            # avoid division by zero
                            maxlen = max(i[3], j[3], 1)
                            minlen = min(i[3], j[3])
                            if float(minlen) / float(maxlen) >= bias:
                                i[-1] = 1
                    except Exception:
                        continue
                else:
                    # DEL/DUP/INV: interval overlap check
                    try:
                        if i[0] != j[0]:
                            continue
                        # i: [chr, pos, END, SVLEN, matched_flag]
                        start_i = i[1]
                        end_i = i[2]
                        start_j = j[1]
                        end_j = j[2]
                        # extend by offect
                        if max(start_i - offect, start_j) <= min(end_i + offect, end_j):
                            # check length ratio (avoid div by zero)
                            maxlen = max(i[3], j[3], 1)
                            minlen = min(i[3], j[3])
                            if float(minlen) / float(maxlen) >= bias:
                                i[-1] = 1
                    except Exception:
                        continue


def statistics_true_positive(callset, SVTYPE):
    """
    Count total records and matched records (matched_flag==1).
    If SVTYPE == "ALL", include all SV types.
    Returns (record, true_record)
    """
    record = 0
    true_record = 0
    if SVTYPE == "ALL":
        for svtype in callset:
            for i in callset[svtype]:
                record += 1
                if i[-1] == 1:
                    true_record += 1
    else:
        if SVTYPE not in callset:
            return record, true_record
        for i in callset[SVTYPE]:
            record += 1
            if i[-1] == 1:
                true_record += 1
    return record, true_record


def main_ctrl(args):
    logging.info("Load SV callset of selected callers.")
    call_child = load_callset(args.F1)
    call_father = load_callset(args.MP)
    call_mother = load_callset(args.FP)

    logging.info("Evaluate matches between child and parents (position/length/breakpoint matching).")

    # Mark parent records that are matched by child (not strictly needed for MDR but retained)
    eva_record(call_child, call_father, args.bias, args.offect)
    eva_record(call_child, call_mother, args.bias, args.offect)

    # Mark child records that are matched by parents (we want to know which child records have matches)
    eva_record(call_father, call_child, args.bias, args.offect)
    eva_record(call_mother, call_child, args.bias, args.offect)

    # Print a clear header with fixed column widths (CONSISTENT column removed)
    header = "{:<6} {:<5} {:>9} {:>12} {:>9}".format("SAMPLE", "SV", "TOTAL", "INCONSISTENT", "MDR(%)")
    logging.info(header)
    logging.info("-" * len(header))

    # Desired output order: DEL, INS, DUP, INV, BND, Total
    svtypes = ["DEL", "INS", "DUP", "INV", "BND", 'ALL']
    for sv in svtypes:
        # Child (F1): total child SVs of this type, and how many were matched by either parent
        record, true_record = statistics_true_positive(call_child, sv)
        if record == 0:
            inconsistent = 0
            mdr = 0.0
        else:
            inconsistent = record - true_record
            mdr = 100.0 * float(inconsistent) / float(record)

        display_sv = 'Total' if sv == 'ALL' else sv
        # Columns: SAMPLE SV TOTAL INCONSISTENT MDR(%) 
        line = "{:<6} {:<5} {:>9d} {:>12d} {:>9.2f}".format('F1', display_sv, record, inconsistent, mdr)
        logging.info(line)


def main(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()
    main_ctrl(args)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))


def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="Trio_eval", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("MP", type=str, help="Male parent callsets (VCF-like)")
    parser.add_argument('FP', type=str, help="Female parent callsets (VCF-like)")
    parser.add_argument('F1', type=str, help="Offspring callsets (VCF-like)")
    parser.add_argument('-b', '--bias', help="Bias of overlapping length similarity [%(default)s]", default=0.7, type=float)
    parser.add_argument('-o', '--offect', help="Offset (bp) allowed for breakpoint/position matching [%(default)s]", default=1000, type=int)
    args = parser.parse_args(argv)
    return args


def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    logging.info("Running %s" % " ".join(sys.argv))


if __name__ == '__main__':
    main(sys.argv[1:])
