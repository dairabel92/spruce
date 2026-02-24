#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Original License:

(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 07 August 2012 21:08 PDT (-0700)

Modifications for SPrUCE:

Edited for SPrUCE 2024 Daira Melendez, Ali Osman Berk Sapci
"""

import os
import re
import sys
import math
import time
import glob
import argparse
import multiprocessing
from collections import Counter, defaultdict
from random import choice
import statistics
import pdb
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.interpolate import UnivariateSpline
import pandas as pd
import logging

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

VERSION="v0.1.0"

def setup_logging(args):
    logging.basicConfig(
        level=getattr(logging, args.verbosity),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    return logging.getLogger("spruce"), "spruce"

def is_dir(path):
    if not os.path.isdir(path):
        raise NotADirectoryError(f"Directory not found: {path}")
    return path

class FullPaths(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(values))

def get_alignment_files(log, directory, file_format):
    valid_ext = ["fasta", "fa", "fas", "nexus", "phy", "phylip", "clustal", "emboss", "stockholm"]
    files = []
    for ext in valid_ext:
        files.extend(glob.glob(os.path.join(directory, f"*.{ext}")))
    if not files:
        log.critical(f"No alignment files found in {directory}.")
        sys.exit(1)
    return files

def get_args():
    parser = argparse.ArgumentParser(
        description=f"SPrUCE ({VERSION}): Sigmoid Pi Requiring Ultraconserved Elements"
    )
    parser.add_argument(
        "--alignments",
        required=True,
        type=is_dir,
        action=FullPaths,
        help="""The directory containing the alignment files""",
    )
    parser.add_argument(
        "--output-file",
        required=True,
        action=FullPaths,
        help="""The name of the CSV file of substitutions""",
    )
    parser.add_argument(
        "--input-format",
        dest="input_format",
        choices=["fasta", "nexus", "phylip", "clustal", "emboss", "stockholm"],
        default="fasta",
        help="""The input alignment format (default: fasta)""",
    )
    parser.add_argument(
        "--cores", type=int, default=1, help="""The number of cores to use (default: 1)""",
    )
    parser.add_argument(
        "--mode",
        type=int,
        default=1,
        help="""1: read alignments; 2: continue, read previous smiley output (default: 1)""",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use (default: INFO)""",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["stack", "concat"],
        default="stack",
        help="""The method to combine UCEs (default: stack)""",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=math.inf,
        help="""The length of the flanking region to consider""",
    )
    parser.add_argument(
        "--min-bases",
        type=int,
        default=0,
        help="""The minimum number of bases (non-gap) for a position to be considered (default: 0)""",
    )
    parser.add_argument(
        "--use-weights",
        default=True,
        type=bool,
        # action=argparse.BooleanOptionalAction,
        help="""Use variance of the esimator to down-weight position with low signal (default: True)""",
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs""",
    )
    return parser.parse_args()


def replace_gaps_at_start_and_ends(seq):
    """walk end from the ends of alignments to locate and replace gaps at ends"""
    begin, end = [], []
    for start, base in enumerate(seq):
        if base == "-":
            begin.append("?")
            continue
        else:
            break
    for stop, base in enumerate(seq[::-1]):
        if base == "-":
            end.append("?")
            continue
        else:
            stop = len(seq) - stop
            break
    newseq = "".join(begin) + str(seq[start:stop]) + "".join(end)
    return Seq(newseq)


def replace_gaps(aln):
    """we need to determine actual starts of alignments"""
    new_aln = MultipleSeqAlignment([])
    for taxon in aln:
        seq = replace_gaps_at_start_and_ends(taxon.seq)
        new_aln.append(
            SeqRecord(seq, id=taxon.id, name=taxon.name, description=taxon.description,)
        )
    return new_aln


def worker(work):
    arguments, f = work
    results = {}
    name_map = {}
    base_count = {}
    locus = os.path.splitext(os.path.basename(f))[0]
    print(f)
    aln = AlignIO.read(f, arguments.input_format)
    # map taxon position in alignment to name
    for idx, taxon in enumerate(aln):
        name_map[idx] = taxon.id
        results[taxon.id] = {
            "insertion": [],
            "deletion": [],
            "substitution": [],
            "majallele": [],
            None: [],
        }
    # get rid of end gappiness, since that makes things a problem
    # for indel ID. Substitute "?" at the 5' and 3' gappy ends.
    # we assume internal gaps are "real" whereas end gaps usually
    # represent missing data.
    aln = replace_gaps(aln)

    positions_with_multiple_alleles = []  # New list to store positions with >2 alleles

    for idx in range(aln.get_alignment_length()):
        col = aln[:, idx]
        col = col.lower()
        # strip the "n" or "N"
        bases = re.sub("N|n|\\?", "", col)
        # count total number of sites considered
        base_count[idx] = len(bases)
        # if all the bases are replace N|n|?, skip
        if len(bases) == 0:
            pass
        # if there is only 1 base, make it the major allele
        elif len(set(bases)) == 1:
            major = bases[0].lower()
        # if there are multiple alleles, pick the major allele
        else:
            # count all the bases in a column
            count = Counter(bases)
            # Check for positions with more than two alleles
            if len(count) > 2:
                positions_with_multiple_alleles.append(idx)
                continue

            # get major allele if possible
            count_of_count = Counter(list(count.values()))
            # we can't have a tie
            if count_of_count[max(count_of_count)] == 1:
                major = count.most_common()[0][0]
            else:
                # randomly select a major allele (excluding gaps)
                # when there is a tie
                common_bases = []
                for base, c in count.most_common():
                    base = base.lower()
                    # bases can be any of IUPAC set except N|n
                    if c == max(count_of_count) and base in [
                        "a",
                        "c",
                        "t",
                        "g",
                        "r",
                        "y",
                        "s",
                        "w",
                        "k",
                        "m",
                        "b",
                        "d",
                        "h",
                        "v",
                    ]:
                        common_bases.append(base)
                # randomly select 1 of the bases
                major = choice(common_bases)
            # now, check for indels/substitutions
        for pos, base in enumerate(col):
            base = base.lower()
            if base in ["N", "n", "?"]:
                results[name_map[pos]][None].append(idx)
            elif base == major:
                results[name_map[pos]]["majallele"].append(idx)
            elif major == "-" and base != "-":
                results[name_map[pos]]["insertion"].append(idx)
            elif base == "-" and major != "-":
                results[name_map[pos]]["deletion"].append(idx)
            elif base != "-" and major != "-":
                results[name_map[pos]]["substitution"].append(idx)
    sys.stdout.write(".")
    sys.stdout.flush()
    return (locus, results, aln.get_alignment_length(), base_count)


def calculate_threshold(uce_data):
    # Calculate ECDF
    Fn = uce_data["bp"].rank(method="first") / len(uce_data)
    out = pd.DataFrame(
        {"bp": np.log2(uce_data["bp"].values), "cdf": np.log2(Fn.values + 1)}
    )

    # Sort the values to ensure they are strictly increasing
    out = out.sort_values(by="bp").reset_index(drop=True)

    # Remove duplicates to ensure strict monotonicity
    out = out.loc[out["bp"].diff() != 0]

    # Smooth the ECDF using cubic spline interpolation
    spline = UnivariateSpline(out["bp"], out["cdf"], s=0, k=3)
    s_bp = np.linspace(min(out["bp"]), max(out["bp"]), num=len(out))
    s_cdf = spline(s_bp)

    # Create the output DataFrame
    out["s_bp"] = s_bp
    out["s_cdf"] = s_cdf

    # Calculate second derivative and find the threshold
    out["s_grad"] = np.gradient(np.gradient(out["s_cdf"]))
    out["c"] = out["s_grad"] < 10 ** -5

    threshold = 2 ** np.max(out["s_bp"][out["s_grad"] < 0])
    return threshold


def estimate_theta(positions, frequencies, ns, bps, args):
    def gompertz2(x, theta, b, c):
        return theta * np.exp(-b * np.exp(-c * np.abs(x)))

    def minimize_function_wstack(parameters):
        x, b, c = parameters
        return sum(
            (n - 1) / (n + 1) * (f - gompertz2(p, x, b, c)) ** 2
            for (p, f, n, bp) in zip(positions, frequencies, ns, bps)
            if bp > args.min_bases
        )

    def minimize_function_wconcat(parameters):
        x, b, c = parameters
        return sum(
            sum(
                (n - 1) / (n + 1) * (f - gompertz2(p, x, b, c)) ** 2
                for f, n, bp in zip(f_list, n_list, bp_list)
                if bp > 1
            )
            for (p, f_list, n_list, bp_list) in zip(positions, frequencies, ns, bps)
        )

    def minimize_function_stack(parameters):
        x, b, c = parameters
        return sum(
            (f - gompertz2(p, x, b, c)) ** 2
            for (p, f, bp) in zip(positions, frequencies, bps)
            if bp > args.min_bases
        )

    def minimize_function_concat(parameters):
        x, b, c = parameters
        return sum(
            sum(
                (f - gompertz2(p, x, b, c)) ** 2
                for f, n, bp in zip(f_list, n_list, bp_list)
                if bp > args.min_bases
            )
            for (p, f_list, n_list, bp_list) in zip(positions, frequencies, ns, bps)
        )

    # max position where bp >= threshold
    if args.method == "stack":
        m = max([p for p, bp in zip(positions, bps) if bp >= args.min_bases], default=1)
    elif args.method == "concat":
        # For concat, just ensure there's at least one base (no thresholding needed)
        m = max([p for p, bp in zip(positions, bps) if len(bps) > 0], default=1)

    # fprime similar to flank
    fprime = min(m, args.flank)  # ensure its bigger than 1
    if args.flank == math.inf:
        args.flank = fprime
    print(f"Calculated fprime: {fprime}")

    x0 = [0.001, 10, 0.009]
    bounds = Bounds([0, 1, 0.004], [0.05, 30, 0.03])

    if args.method == "stack":
        if args.use_weights:
            return minimize(
                minimize_function_wstack, x0, bounds=bounds, method="trust-constr"
            )
        else:
            return minimize(
                minimize_function_stack, x0, bounds=bounds, method="trust-constr"
            )
    elif args.method == "concat":
        if args.use_weights:
            return minimize(
                minimize_function_wconcat, x0, bounds=bounds, method="trust-constr"
            )
        else:
            return minimize(
                minimize_function_concat, x0, bounds=bounds, method="trust-constr"
            )
    else:
        raise ValueError("Only 'stack' and 'concat' are available for 'method'.")


def main():
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    if args.mode == 1:
        # create database name
        # iterate through all the files to determine the longest alignment
        log.info("Checking alignments....")
        work = [
            (args, f)
            for f in get_alignment_files(log, args.alignments, args.input_format)
        ]
        if args.cores > 1:
            pool = multiprocessing.Pool(args.cores)
            results = pool.map(worker, work)
        else:
            results = list(map(worker, work))
        print("")
        # fill the individual/locus/position specific table
        if args.method == "stack":
            rsf = dict()
            rsn = dict()
            cbp = dict()
        elif args.method == "concat":
            rsf = defaultdict(list)
            rsn = defaultdict(list)
            cbp = defaultdict(list)
            ruce = defaultdict(list)
        else:
            raise ValueError("Only 'stack' and 'concat' are available for 'method'.")
        rcf = dict()

        centers = {}
        # adding locus set for thresold
        locus_set = set()
        for locus, result, length, bases in results:
            # get approximate center of alignment
            center = centers.get(locus, length // 2)
            if abs(center - length // 2) > 50:
                continue
            # fill locus table
            # we also want a locus specific list of all variable positions
            # NOTE: currently only doing this for substitutions
            maj, subs, dels, ins, n = [], [], [], [], []
            # get all substitution locations across individuals
            for k, v in list(result.items()):
                subs.extend(v["substitution"])
                dels.extend(v["deletion"])
                ins.extend(v["insertion"])
                n.extend(v[None])
            # get a count of variability by position in BP
            subs_cnt = Counter(subs)
            dels_cnt = Counter(dels)
            ins_cnt = Counter(ins)
            n_cnt = Counter(n)
            # iterate over counts of all positions - having subs and not having subs
            # then add those + any sub location to the DB
            for pos in sorted(bases.keys()):
                substitutions = subs_cnt[pos]
                bsp = bases[pos]
                insertions = ins_cnt[pos]
                deletions = dels_cnt[pos]
                missing = n_cnt[pos]
                poscen = pos - center
                n = bsp - insertions - deletions
                if ((n) > 1) and (abs(poscen) < args.flank):
                    if args.method == "stack":
                        fij = (
                            2.0
                            * (substitutions)
                            * (n - substitutions)
                            / (n)
                            / (n - 1.0)
                        )
                        if fij < 0:
                            raise ValueError( "This should not happen %d, %d, %d, %d, %d, %d, %d, %d" %(n, substitutions, bsp , insertions , deletions , missing, poscen, locus))

                        rsf[poscen] = 0.0 + rsf.get(poscen, 0) + fij
                        cbp[poscen] = cbp.get(poscen, 0) + bsp
                        rsn[poscen] = rsn.get(poscen, 0) + n
                    elif args.method == "concat":
                        fij = (
                            2 * (substitutions) * (n - substitutions) / (n) / (n - 1.0)
                        )
                        rsf[poscen].append(fij)
                        cbp[poscen].append(bsp)
                        rsn[poscen].append(n)
                        ruce[poscen].append(locus)
                    else:
                        raise ValueError(
                            "Only 'stack' and 'concat' are available for 'method'."
                        )
                    locus_set.add(locus)
                    rcf[poscen] = rcf.get(poscen, 0) + 1

        if args.method == "stack":
            # Take averages
            for k, v in rsf.items():
                rsf[k] = rsf[k] / rcf[k]
            outf = open(args.output_file, "w")
            outf.write("distance_from_center,freq,bp,n\n")
            for k in sorted(rsf.keys()):
                outf.write("%s,%s,%s,%s\n" % (k, rsf[k], cbp[k], rsn[k]))
            outf.close()
        elif args.method == "concat":
            outf = open(args.output_file, "w")
            outf.write("distance_from_center,freq,bp,n\n")
            for k in rsf.keys():
                for i in range(len(rsf[k])):
                    outf.write(
                        "%s,%s,%s,%s,%s\n"
                        % (ruce[k][i], k, rsf[k][i], cbp[k][i], rsn[k][i])
                    )
            outf.close()
        else:
            raise ValueError("Only 'stack' and 'concat' are available for 'method'.")
    else:
        if args.method == "stack":
            rsf = dict()
            cbp = dict()
            rsn = dict()
            for line in open(args.output_file).readlines()[1:]:
                r = line.strip().split(",")
                rsf[r[0]] = float(r[1])
                cbp[r[0]] = float(r[2])
                rsn[r[0]] = float(r[3])
        elif args.method == "concat":
            rsf = defaultdict(list)
            cbp = defaultdict(list)
            rsn = defaultdict(list)
            for line in open(args.output_file).readlines()[1:]:
                r = line.strip().split(",")
                rsf[r[0]].append(float(r[1]))
                cbp[r[0]].append(float(r[2]))
                rsn[r[0]].append(float(r[3]))
        else:
            raise ValueError("Only 'stack' and 'concat' are available for 'method'.")

    if args.method == "stack" and args.min_bases == 0:
        uce_data_df = pd.DataFrame({"bp": cbp.values()})
        threshold = calculate_threshold(uce_data_df)
        if abs(threshold - uce_data_df["bp"].mean()) > (
            uce_data_df["bp"].std() * 2
        ) or (threshold > uce_data_df["bp"].mean()):
            threshold = len(locus_set) * 2
        print("Calculated threshold for min-base: {}".format(threshold))
        args.min_bases = threshold
    elif args.method == "concat":
        args.min_base = 1

    start_time = time.time()

    res = estimate_theta(
        rsf.keys(), rsf.values(), rsn.values(), [cbp[k] for k in rsf.keys()], args
    )

    end_time = time.time()
    elapsed_time = end_time - start_time

    print("The estimated theta is: " + str(res.x[0]))
    print("All Gompertz parameters are: " + str(res.x))
    print(f"Output time: {elapsed_time:.2f} seconds")  # yay print out elapsed time


if __name__ == "__main__":
    main()
