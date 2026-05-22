#!/usr/bin/env python

import argparse
import numpy as np
from fasta import Sequence

STOP, DIAG, UP, LEFT = 0, 1, 2, 3
M_SOURCE, X_SOURCE, Y_SOURCE = 0, 1, 2

def scoring_matrix(seq1: str, seq2: str, pam=None, blosum=None, mode="global", match=2, mismatch=-1, go=-2, ge=-1):
    """Claculates the scorint matrix for two candidate sequences for pairwise alignment
    :param seq1: str: string of nucleotide sequences for the first DNA sequence
    :param seq2: str: string of nucleotide sequences for the second DNA sequence
    :param mode: str: alignment node. options are local and global
    :param match int: score for a matching base between the two sequences, default is 1
    :param mismatch int: penalty for taking a mismatch to the next cell, default is -1
    :param go int: penalty for opening a gap, default is -3
    :param ge int: penealty for gap extension, default is -2
    """
    if mode not in {"global","local"}:
        raise ValueError("only valid entries for mode are global and local")

    n, m = len(seq1), len(seq2)

    M=np.full((n+1, m+1), -np.inf)
    X=np.full((n+1, m+1), -np.inf)
    Y=np.full((n+1, m+1), -np.inf)

    trace_M = np.full((n+1, m+1), -1, dtype=int)
    trace_X = np.full((n+1, m+1), -1, dtype=int)
    trace_Y = np.full((n+1, m+1), -1, dtype=int)


    #global
    if mode == "global":
        M[0,0] = 0
        for i in range(1, n+1):
            X[i,0]=go + i * ge
            trace_X[i,0] = X_SOURCE
        for j in range(1, m+1):
            Y[0,j]=go + j * ge
            trace_Y[0,j] = Y_SOURCE

    for i in range(1, n+1):
        for j in range(1, m+1):
            if pam is not None and blosum is not None:
                raise ValueError("Use either PAM or BLOSUM, not both")

            if pam is not None:
                s = pam[seq1[i-1], seq2[j-1]]
            elif blosum is not None:
                s = blosum[seq1[i-1], seq2[j-1]]
            else:
                s = match if seq1[i-1] == seq2[j-1] else mismatch

            if mode == "global":
                candidates_M = [M[i-1,j-1], X[i-1,j-1], Y[i-1,j-1]]
                best_prev_M = max(candidates_M)
                M[i,j]= best_prev_M + s
                trace_M[i,j] = candidates_M.index(best_prev_M)

                candidates_X = [M[i-1,j] + go + ge, X[i-1,j] + ge]
                best_prev_X = max(candidates_X)
                X[i,j]=best_prev_X
                trace_X[i,j] = candidates_X.index(best_prev_X)

                candidates_Y = [M[i,j-1] + go + ge, Y[i,j-1] + ge]
                best_prev_Y = max(candidates_Y)
                Y[i,j]=best_prev_Y
                trace_Y[i,j] = candidates_Y.index(best_prev_Y)

            elif mode == "local":
                candidates_M = [0, M[i-1,j-1] + s, X[i-1,j-1] + s, Y[i-1,j-1] + s]
                M[i,j]= max(candidates_M)
                trace_M[i,j] = candidates_M.index(M[i,j])

                candidates_X = [0, M[i-1,j] + go + ge, X[i-1,j] + ge]
                X[i,j]= max(candidates_X)
                trace_X[i,j] = candidates_X.index(X[i,j])


                candidates_Y = [0, M[i,j-1] + go + ge, Y[i,j-1] + ge]
                Y[i,j]= max(candidates_Y)
                trace_Y[i,j] = candidates_Y.index(Y[i,j])

    return M, X, Y, trace_M, trace_X, trace_Y

def back_track(seq1, seq2, M, X, Y, trace_M, trace_X, trace_Y, mode="global"):
    aln1 = []
    aln2 = []

    n, m = len(seq1), len(seq2)

    if mode == "global":
        end_scores = [M[n,m],X[n,m],Y[n,m]]
        state= end_scores.index(max(end_scores))
        i, j = n, m

    elif mode == "local":
        max_score = -np.inf
        state = None

        for mat, s in zip([M, X, Y], [M_SOURCE, X_SOURCE, Y_SOURCE]):
            idx = np.unravel_index(np.argmax(mat), mat.shape)
            if mat[idx] > max_score:
                max_score = mat[idx]
                i, j = idx
                state = s

    while True:
        if mode == "global" and i == 0 and j ==0:
            break
        if mode == "local":
            if state == M_SOURCE and M[i,j] == 0:
                break
            if state == X_SOURCE and X[i,j] == 0:
                break
            if state == Y_SOURCE and Y[i,j] == 0:
                break

        if state == M_SOURCE:
            #came from diagonal
            prev_state = trace_M[i,j]
            # state emission is a match or a mismatch
            aln1.append(seq1[i-1])
            aln2.append(seq2[j-1])
            i -= 1
            j -= 1
            state = prev_state
        elif state == X_SOURCE:
            #gap in seq2 (vertical move)
            prev_state = trace_X[i,j]
            aln1.append(seq1[i-1])
            aln2.append("-")
            i -= 1
            state = prev_state
        elif state == Y_SOURCE:
            #gap in seq1 (horizontal move)
            prev_state = trace_Y[i,j]
            aln1.append("-")
            aln2.append(seq2[j-1])
            j -= 1
            state = prev_state
    return "".join(reversed(aln1)), "".join(reversed(aln2))

def use_sub_mat(path):
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    header = lines[0].split()
    symbols = header
    scores = {}

    for line in lines[1:]:
        parts = line.split()
        row_symbol = parts[0]
        values = list(map(int, parts[1:]))

        for col_symbol, val in zip(symbols, values):
            scores[(row_symbol, col_symbol)] = val

    return scores


def parse_args():
    parser = argparse.ArgumentParser(description="Sequence alignment tool for DNA sequences")
    parser.add_argument("ifile",help="Input file containing DNA/Protein sequences in fasta format(.fa, .fas, .fasta)")

    subparse = parser.add_subparsers(dest="mode",required=True)

    # -----------Pairwise alignment---------
    pairwise = subparse.add_parser("pairwise",help="Pairwise Sequence alignment")
    pairwise.add_argument("--needle",action="store_true",help="Global alignment (Needleman-Wunsch)")
    pairwise.add_argument("--smith",action="store_true",help="Local alignment (Smith-Waterman)")
    pairwise.add_argument("--seq1",required=True,help="ID of first sequence")
    pairwise.add_argument("--seq2",required=True,help="ID of second sequence")
    pairwise.add_argument("--dna", help="DNA alignment using PAM substitution matrix")
    pairwise.add_argument("--prot", help="Protein alignment using BLOSUM substitution matrix")

    # ----------Multiple sequence alignment (MSA)-------
    msa = subparse.add_parser("msa",help="Multiple sequence alignment")
    msa.add_argument("--pam", help="Path to PAM matrix")
    msa.add_argument("--blosum", help="Path to BLOSUM matrix")

    return parser.parse_args()

def main():
    args = parse_args()
    sequences = Sequence.fasta_file(args.ifile)
    print("Initiated")
    seq1 = sequences[args.seq1].seq
    seq2 = sequences[args.seq2].seq
    if args.needle:
        mode = "global"
    elif args.smith:
        mode = "local"
    pam = None
    blosum = None
    if args.dna:
        pam = use_sub_mat(args.dna)

    elif args.prot:
        blosum= use_sub_mat(args.prot)

    M, X, Y, trace_M, trace_X, trace_Y = scoring_matrix(seq1, seq2, pam=pam, blosum=blosum,mode=mode)
    aln1, aln2 = back_track(seq1, seq2, M, X, Y, trace_M, trace_X, trace_Y, mode=mode)
    print(aln1)
    print(aln2)
    print("Finished!")

if __name__ == "__main__":
    main()
