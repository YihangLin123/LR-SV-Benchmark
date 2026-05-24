#!/usr/bin/env python3
import sys
import argparse
import logging
import time
import re

def norm_chr(chrom):
    if chrom is None:
        return chrom
    chrom = str(chrom)
    if chrom.lower().startswith('chr'):
        return chrom[3:]
    return chrom

def pase_info(seq):
    info = {'SVLEN': 0, 'END': 0, 'POS2': 0, "SVTYPE": '', "RE": 0, "CHR2": ''}
    for item in seq.split(';'):
        if '=' not in item:
            continue
        k, v = item.split('=', 1)
        if k in ["SVLEN", "END", "RE", "POS2"]:
            try:
                info[k] = abs(int(float(v)))
            except:
                pass
        elif k == "CHR2":
            info[k] = v
        elif k == "SVTYPE":
            info[k] = v  
    return info

def parse_bnd_alt(alt):
    if not alt:
        return None, None, None
    m = re.search(r'([A-Za-z0-9_.]+):(\d+)', alt)
    if not m:
        return None, None, None
    chr2 = m.group(1)
    pos2 = int(m.group(2))

    if '[' in alt:
        form = '[[N' if alt.startswith('[') else 'N[['
    elif ']' in alt:
        form = ']]N' if alt.startswith(']') else 'N]]'
    else:
        form = 'N[['
    return chr2, pos2, form

def safe_chr_cmp(a, b):
    try:
        return int(a) <= int(b)
    except:
        return str(a) <= str(b)

def phase_GT(sample_field, format_field):
    try:
        fmt = format_field.split(':')
        sample_vals = sample_field.split(':')
        gt_index = fmt.index('GT')
        gt = sample_vals[gt_index]
    except Exception:
        gt = sample_field.split(':')[0] if sample_field else './.'
    if gt in ['0/1', '1/0', '0|1', '1|0']:
        return 'het'
    elif gt in ['1/1', '1|1']:
        return 'hom'
    else:
        return 'unknown'

def load_callset(path, svtype_list):
    callset = dict()
    abtype = dict()
    if not path:
        return callset, abtype
    with open(path, 'r') as file:
        for line in file:
            seq = line.strip('\n').split('\t')
            if len(seq) < 5 or seq[0].startswith('#'):
                continue

            chr_raw = seq[0]
            try:
                pos = int(seq[1])
            except:
                continue
            info = pase_info(seq[7])

            if len(svtype_list) == 3 and info['SVTYPE'] == "DUP":
                info['SVTYPE'] = "INS"

            if info['SVTYPE'] in svtype_list:
                chr_norm = norm_chr(chr_raw)
                format_field = seq[8] if len(seq) > 8 else ''
                sample_field = seq[9] if len(seq) > 9 else ''

                if info['SVTYPE'] == "BND":
                    chr2_alt, pos2_alt, form = parse_bnd_alt(seq[4])
                    if not chr2_alt: continue
                    chr2_norm = norm_chr(chr2_alt)
                    
                    if info['SVTYPE'] not in callset:
                        callset[info['SVTYPE']] = list()

                    
                    if safe_chr_cmp(chr_norm, chr2_norm):
                        callset[info['SVTYPE']].append([chr_norm, pos, chr2_norm, pos2_alt, form, phase_GT(sample_field, format_field), 0])
                    else:
                        callset[info['SVTYPE']].append([chr2_norm, pos2_alt, chr_norm, pos, form, phase_GT(sample_field, format_field), 0])
                else:
                    if info['SVTYPE'] not in callset:
                        callset[info['SVTYPE']] = list()
                    slen = info['SVLEN'] if info['SVLEN'] != 0 else abs(info.get('END', pos) - pos)
                    callset[info['SVTYPE']].append([norm_chr(chr_raw), pos, info.get('END', 0), slen, phase_GT(sample_field, format_field), 0])
            else:
                tkey = info.get('SVTYPE', 'UNK')
                abtype[tkey] = abtype.get(tkey, 0) + 1
    return callset, abtype

def eval(call, ans, bias, offect, opt, genotype):
    for svtype in call:
        if svtype not in ans:
            continue
        for i in call[svtype]:
            for j in ans[svtype]:
                if i[0] != j[0]: continue
                
                if svtype == 'INS':
                    if abs(i[1] - j[1]) <= offect and float(min(i[3],j[2])/max(i[3],j[2])) >= bias:
                        i[-1] = 1; j[2+opt] = 1
                        if i[4] == genotype.get(j[0], ''):
                            i[-1] = 2; j[2+opt] = 2
                            
                elif svtype == 'BND':
                   
                    if i[2] == j[2] and i[4] == j[4]: 
                        if abs(i[1]-j[1]) <= offect and abs(i[3]-j[3]) <= offect:
                            i[-1] = 1; j[4+opt] = 1
                            if i[5] == genotype.get(j[0], '') or i[5] == genotype.get(j[2], ''):
                                i[-1] = 2; j[4+opt] = 2
                else:
                    if max(i[1]-offect, j[1]) <= min(i[2]+offect, j[2]) and float(min(i[3],j[3])/max(i[3],j[3])) >= bias:
                        i[-1] = 1; j[3+opt] = 1
                        if i[4] == genotype.get(j[0], ''):
                            i[-1] = 2; j[3+opt] = 2

def statistics(call, ans, opt, res):
    for svtype in call:
        tp = sum(1 for e in call[svtype] if e[-1] >= res)
        logging.info('TP-%d of %s:\t%d\t%d' % (res, svtype, tp, len(call[svtype])))
    for svtype in ans:
        idx = 2+opt if svtype=='INS' else (4+opt if svtype=='BND' else 3+opt)
        fn = sum(1 for e in ans[svtype] if e[idx] >= res)
        logging.info('TN-%d of %s:\t%d\t%d' % (res, svtype, fn, len(ans[svtype])))

typetrans = {'insertion':'INS', 'deletion':'DEL', 'inversion':'INV', 'tandem duplication':'DUP', 'reciprocal translocation':'BND'}

def load_ans(path, ntools):
    ansbed = dict()
    min_cols = max(12, 6 + ntools) 
    with open(path, 'r') as file:
        for line in file:
            seq = line.strip('\n').split('\t')
            if len(seq) < 4: continue
            chr1 = norm_chr(seq[0])
            svtype = typetrans.get(seq[3], 'BND')
            start, end = int(seq[1]), int(seq[2])
            if svtype not in ansbed: ansbed[svtype] = []
            
            if svtype == 'BND':
                parts = seq[4].split(':')
                chr2, start2 = norm_chr(parts[1]), int(parts[2])
                s1, s2 = parts[3][0], parts[4][0]
                length = end - start
                
                conns = []
                if s1 == 'f' and s2 == 'f':
                    conns = [("N[[", chr1, start, chr2, start2), ("]]N", chr1, start, chr2, start2), 
                             ("]]N", chr1, end, chr2, start2+length), ("N[[", chr1, end, chr2, start2+length)]
                elif s1 == 'f' and s2 == 'r':
                    conns = [("N[[", chr1, start, chr2, start2), ("[[N", chr1, start, chr2, start2+length), 
                             ("N]]", chr1, end, chr2, start2), ("]]N", chr1, end, chr2, start2+length)]
                elif s1 == 'r' and s2 == 'f':
                    conns = [("N]]", chr1, start, chr2, start2+length), ("]]N", chr1, start, chr2, start2), 
                             ("[[N", chr1, end, chr2, start2), ("N[[", chr1, end, chr2, start2+length)]
                else: # r and r
                    conns = [("N]]", chr1, start, chr2, start2+length), ("N]]", chr1, end, chr2, start2), 
                             ("[[N", chr1, end, chr2, start2), ("[[N", chr1, start, chr2, start2+length)]
                
                for f, c1, p1, c2, p2 in conns:
                    if safe_chr_cmp(c1, c2):
                        ansbed[svtype].append([c1, p1, c2, p2, f] + [0]*min_cols)
                    else:
                        ansbed[svtype].append([c2, p2, c1, p1, f] + [0]*min_cols)
            elif svtype == 'INS':
                ansbed[svtype].append([chr1, start, len(seq[4]) if len(seq)>4 else 0] + [0]*min_cols)
            else:
                ansbed[svtype].append([chr1, start, end, end-start+1] + [0]*min_cols)
    return ansbed

def load_gt(path):
    GT = dict()
    with open(path) as fh:
        for line in fh:
            seq = line.strip().split('\t')
            if len(seq) < 2: continue
            try:
                score = float(seq[-1])
            except: score = 0.0
            GT[norm_chr(seq[0])] = 'hom' if score > 80 else ('het' if score > 20 else 'None')
    return GT

def main_ctrl(args):
    ntools = len(args.tools) if args.tools else 0
    ans = load_ans(args.ans, ntools)
    genotype = load_gt(args.gt)
    
    svtypes = {"BND":["BND"], "DUP":["INS","DUP"], "IID":["INS","INV","DEL"]}.get(args.choice)

    for idx, path in enumerate(args.tools or []):
        opt_label = idx + 1
        call, ab = load_callset(path, svtypes)
        logging.info(f"Evaluation on {path} callsets...")
        eval(call, ans, args.bias, args.offect, opt_label, genotype)
        statistics(call, ans, opt_label, 1)
        statistics(call, ans, opt_label, 2)

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("choice", help="IID/DUP/BND")
    parser.add_argument("ans", help="Ground truth BED")
    parser.add_argument("gt", help="Genotype file")
    parser.add_argument('-t', '--tools', nargs='+', help="VCF paths")
    parser.add_argument('-b', '--bias', default=0.7, type=float)
    parser.add_argument('-o', '--offect', default=1000, type=int)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    main_ctrl(args)

if __name__ == '__main__':
    main(sys.argv[1:])