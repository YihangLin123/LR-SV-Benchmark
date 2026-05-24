#!/usr/bin/env python3
import sys
import argparse
import logging
import time
import re

def norm_chr(chrom):
    """统一染色体命名：去掉 'chr' 前缀并返回（保持 None 安全）。"""
    if chrom is None:
        return chrom
    chrom = str(chrom)
    if chrom.lower().startswith('chr'):
        return chrom[3:]
    return chrom

def pase_info(seq):
    """从 INFO 字段提取常用键；支持 POS2（若有）、CHR2、END、SVLEN、SVTYPE、RE。"""
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
            info[k] = v  # 取完整值，不截断
    return info

def parse_bnd_alt(alt):
    """
    从 ALT 字段（如 N]3:76470743] 或 T]3:76470743] 等）提取 chr2, pos2, form。
    如果解析失败，返回 (None, None, None)。
    """
    if not alt:
        return None, None, None
    m = re.search(r'([A-Za-z0-9_.]+):(\d+)', alt)
    if not m:
        return None, None, None
    chr2 = m.group(1)
    pos2 = int(m.group(2))

    # 简单判断 form：根据方括号在 ALT 中相对位置
    if '[' in alt and ']' not in alt:
        if alt.find('[') < alt.find(chr2):
            form = '[[N'
        else:
            form = 'N[['
    elif ']' in alt and '[' not in alt:
        if alt.find(']') < alt.find(chr2):
            form = ']]N'
        else:
            form = 'N]]'
    else:
        if alt.find('[') != -1:
            form = 'N[['
        elif alt.find(']') != -1:
            form = 'N]]'
        else:
            form = 'N[['
    return chr2, pos2, form

def safe_chr_cmp(a, b):
    """尽量稳定地比较两个染色体名的顺序（用于决定哪端先写）。"""
    try:
        ai = int(a)
        bi = int(b)
        return ai <= bi
    except:
        return str(a) <= str(b)

def phase_GT(sample_field, format_field):
    """更稳健地从 FORMAT 和样本字段里提取 GT 并返回 'het'/'hom'/'unknown'"""
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
            if len(seq) < 5:
                continue
            if seq[0].startswith('#'):
                continue

            chr_raw = seq[0]
            try:
                pos = int(seq[1])
            except:
                continue
            info = pase_info(seq[7])

            # 如果用户选择了 DUP 而调用方报告 DUP，转成 INS（保留你原来的逻辑）
            if len(svtype_list) == 3 and info['SVTYPE'] == "DUP":
                info['SVTYPE'] = "INS"

            if info['SVTYPE'] in svtype_list:
                chr_norm = norm_chr(chr_raw)
                format_field = seq[8] if len(seq) > 8 else ''
                sample_field = seq[9] if len(seq) > 9 else ''

                if info['SVTYPE'] == "BND":
                    # 优先使用 INFO 中的 CHR2/POS2（如果存在），否则从 ALT 解析
                    chr2_info = info.get('CHR2') or None
                    pos2_info = info.get('POS2') or None
                    form = None

                    if not chr2_info or not pos2_info:
                        chr2_alt, pos2_alt, form_alt = parse_bnd_alt(seq[4])
                        if chr2_alt:
                            chr2_info = chr2_info or chr2_alt
                            pos2_info = pos2_info or pos2_alt
                            form = form_alt

                    # 有些工具可能把对端放在 END 字段
                    if not pos2_info and info.get('END', 0) > 0:
                        pos2_info = info.get('END')

                    if not chr2_info or not pos2_info:
                        logging.debug(f"Skipping BND (cannot find CHR2/POS2): {seq[0]}:{seq[1]} ALT={seq[4]} INFO={seq[7]}")
                        tkey = info.get('SVTYPE', 'BND')
                        if tkey not in abtype:
                            abtype[tkey] = 0
                        abtype[tkey] += 1
                        continue

                    chr2_norm = norm_chr(chr2_info)
                    if form is None:
                        _, _, fm = parse_bnd_alt(seq[4])
                        form = fm or 'N[['

                    if info['SVTYPE'] not in callset:
                        callset[info['SVTYPE']] = list()

                    try:
                        if safe_chr_cmp(chr_norm, chr2_norm):
                            callset[info['SVTYPE']].append([chr_norm, pos, chr2_norm, pos2_info, form, phase_GT(sample_field, format_field), 0])
                        else:
                            callset[info['SVTYPE']].append([chr2_norm, pos2_info, chr_norm, pos, form, phase_GT(sample_field, format_field), 0])
                    except Exception as e:
                        logging.debug(f"Fallback append BND due to {e}")
                        callset[info['SVTYPE']].append([chr_norm, pos, chr2_norm, pos2_info, form, phase_GT(sample_field, format_field), 0])

                else:
                    if info['SVTYPE'] not in callset:
                        callset[info['SVTYPE']] = list()
                    if info['SVLEN'] == 0:
                        if info.get('END', 0) > 0:
                            info['SVLEN'] = info['END'] - pos + 1
                        else:
                            info['SVLEN'] = 0
                    callset[info['SVTYPE']].append([norm_chr(chr_raw), pos, info.get('END', 0), info['SVLEN'], phase_GT(sample_field, format_field), 0])
            else:
                tkey = info.get('SVTYPE', 'UNK')
                if tkey not in abtype:
                    abtype[tkey] = 0
                abtype[tkey] += 1
    return callset, abtype

def eval(call, ans, bias, offect, opt, genotype):
    for svtype in call:
        if svtype not in ans:
            if svtype == 'INS':
                for i in call[svtype]:
                    for key in ans:
                        for j in ans[key]:
                            if i[0] == j[0]:
                                if abs(i[1] - j[1]) <= offect and float(min(i[3], j[3]) / max(i[3], j[3])) >= bias:
                                    i[-1] = 1
                                    j[3+opt] = 1
                                    if i[4] == genotype.get(j[0], ''):
                                        i[-1] = 2
                                        j[3+opt] = 2
        else:
            for i in call[svtype]:
                for j in ans[svtype]:
                    if i[0] != j[0]:
                        continue
                    else:
                        if svtype in ['INS']:
                            if abs(i[1] - j[1]) <= offect and float(min(i[3], j[2]) / max(i[3], j[2])) >= bias:
                                j[2+opt] = 1
                                i[-1] = 1
                                if i[4] == genotype.get(j[0], ''):
                                    j[2+opt] = 2
                                    i[-1] = 2
                        elif svtype == 'BND':
                            # BND: i = [chr1,pos1,chr2,pos2,form,gt,flag]
                            if i[2] != j[2]:
                                continue
                            else:
                                if abs(i[1] - j[1]) <= offect and abs(i[3] - j[3]) <= offect:
                                    i[-1] = 1
                                    j[4+opt] = 1
                                    if i[5] == genotype.get(j[0], '') or i[5] == genotype.get(j[2], ''):
                                        i[-1] = 2
                                        j[4+opt] = 2
                        else:
                            if max(i[1]-offect, j[1]) <= min(i[2]+offect, j[2]) and float(min(i[3], j[3]) / max(i[3], j[3])) >= bias:
                                j[3+opt] = 1
                                i[-1] = 1
                                if i[4] == genotype.get(j[0], ''):
                                    j[3+opt] = 2
                                    i[-1] = 2

def statistics(call, ans, opt, res):
    for svtype in call:
        tp = 0
        total = 0
        for ele in call[svtype]:
            total += 1
            if ele[-1] >= res:
                tp += 1
        logging.info('TP-%d of %s:\t%d\t%d' % (res, svtype, tp, total))

    for svtype in ans:
        fn = 0
        total = 0
        for ele in ans[svtype]:
            total += 1
            if svtype == 'INS':
                if ele[2+opt] >= res:
                    fn += 1
            elif svtype == 'BND':
                if ele[4+opt] >= res:
                    fn += 1
            else:
                if ele[3+opt] >= res:
                    fn += 1
        logging.info('TN-%d of %s:\t%d\t%d' % (res, svtype, fn, total))

typetrans = {'insertion':'INS',
             'deletion':'DEL',
             'inversion':'INV',
             'tandem duplication':'DUP',
             'reciprocal translocation':'BND'
             }

def load_ans(path, ntools=0):
    """
    读取 ground-truth bed（或自定义格式）。
    为每条记录在末尾预留一定数量的 0 列以便写入不同工具的评估结果。
    ntools 用于决定预留列数。
    """
    # 确定每行应当预留的列数：确保 indices 2+opt / 3+opt / 4+opt 都不会越界
    # opt 最大为 ntools，因此需要至少 (4 + ntools + 1) 长度，取和原先 12 的最大值
    min_cols = max(12, 5 + ntools)
    ansbed = dict()
    with open(path) as fh:
        for line in fh:
            seq = line.strip().split('\t')
            if len(seq) < 4:
                continue
            chr_raw = seq[0]
            svtype_name = seq[3]
            if svtype_name not in typetrans:
                continue
            svtype = typetrans[svtype_name]
            try:
                start = int(seq[1])
                end = int(seq[2])
            except:
                continue

            if svtype not in ansbed:
                ansbed[svtype] = list()

            if svtype == 'INS':
                # 确保 seq[4] 存在；否则把长度设成 0
                inslen = len(seq[4]) if len(seq) > 4 else 0
                ansbed[svtype].append([norm_chr(chr_raw), start, inslen] + [0]*min_cols)
            elif svtype == 'BND':
                # 期望 seq[4] 中包含 chr:pos:... 的结构（原脚本的假设）
                parts = seq[4].split(':') if len(seq) > 4 else []
                if len(parts) >= 3:
                    chr2 = parts[1]
                    start2 = int(parts[2])
                    strand1 = parts[3] if len(parts) > 3 else 'f'
                    strand2 = parts[4] if len(parts) > 4 else 'f'
                else:
                    # 如果格式不对，跳过该条
                    continue

                chr1n = norm_chr(chr_raw)
                chr2n = norm_chr(chr2)
                length = end - start

                # 根据 strand 情况扩展多个表示（与原脚本行为一致）
                if strand1[0] == 'f':
                    if strand2[0] == 'f':
                        ansbed[svtype].append([chr1n, start, chr2n, start2, "N[["] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, end, chr2n, start2+length, "N[["] + [0]*min_cols)
                    else:
                        ansbed[svtype].append([chr1n, start, chr2n, start2, "N[["] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, start, chr2n, start2+length, "[[N"] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, end, chr2n, start2, "N]]"] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, end, chr2n, start2+length, "]]N"] + [0]*min_cols)
                else:
                    if strand2[0] == 'f':
                        ansbed[svtype].append([chr1n, start, chr2n, start2+length, "N]]"] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, start, chr2n, start2, "]]N"] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, end, chr2n, start2, "[[N"] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, end, chr2n, start2+length, "N[["] + [0]*min_cols)
                    else:
                        ansbed[svtype].append([chr1n, start, chr2n, start2+length, "N]]"] + [0]*min_cols)
                        ansbed[svtype].append([chr1n, end, chr2n, start2, "N]]"] + [0]*min_cols)
            else:
                ansbed[svtype].append([norm_chr(chr_raw), start, end, end-start+1] + [0]*min_cols)  # DEL/INV/DUP
    return ansbed

def load_gt(path):
    GT = dict()
    with open(path) as fh:
        for line in fh:
            seq = line.strip().split('\t')
            if len(seq) < 2:
                continue
            chr_raw = seq[0]
            try:
                score = float(seq[-1])
            except:
                score = 0.0
            if score > 80.0:
                GT[norm_chr(chr_raw)] = 'hom'
            elif 80.0 >= score > 20.0:
                GT[norm_chr(chr_raw)] = 'het'
            else:
                GT[norm_chr(chr_raw)] = 'None'
    return GT

def main_ctrl(args):
    # args.tools may be None or a list of paths
    ntools = len(args.tools) if args.tools else 0

    # load answers with enough columns for all tools
    ans = load_ans(args.ans, ntools)
    genotype = load_gt(args.gt)

    def load_and_eval(path, svtypes, opt_label):
        if path:
            call, ab = load_callset(path, svtypes)
            logging.info(f"The number of calls within abnormal SV type in {path}:")
            for key in ab:
                logging.info(f"<{path}-{key}>\t{ab[key]}.")
            logging.info(f"Evaluation on {path} callsets...")
            eval(call, ans, args.bias, args.offect, opt_label, genotype)
            statistics(call, ans, opt_label, 1)
            statistics(call, ans, opt_label, 2)

    if args.choice == "BND":
        svtypes = ["BND"]
    elif args.choice == "DUP":
        svtypes = ["INS", "DUP"]
    elif args.choice == "IID":
        svtypes = ["INS", "INV", "DEL"]
    else:
        logging.error(f"Unknown choice: {args.choice}")
        return

    if args.tools:
        logging.info("Tool index mapping (opt_label -> tool path):")
        for idx, path in enumerate(args.tools):
            logging.info(f"{idx+1}\t{path}")
        for idx, path in enumerate(args.tools):
            opt_label = idx + 1
            logging.info(f"Processing tool #{opt_label}: {path}")
            load_and_eval(path, svtypes, opt_label)
    else:
        logging.warning("No --tools provided. Nothing to evaluate.")

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="Trio_eval", description="", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("choice", type=str, help="Chose specific SV type.[IID/DUP/BND]")
    parser.add_argument("ans", type=str, help="Ground truth of simulations.")
    parser.add_argument("gt", type=str, help="Genotype in each chromosome.")
    parser.add_argument('-t', '--tools', nargs='+', help="Paths to tool VCFs (order matters). Provide one or more VCFs separated by spaces.", default=None)
    parser.add_argument('-b', '--bias', help="Bias of overlaping.[%(default)s]", default=0.7, type=float)
    parser.add_argument('-o', '--offect', help="Offect of translocation overlaping.[%(default)s]", default=1000, type=int)
    args = parser.parse_args(argv)
    return args

def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    logging.info("Running %s" % " ".join(sys.argv))

def main(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()
    main_ctrl(args)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))

if __name__ == '__main__':
    main(sys.argv[1:])
