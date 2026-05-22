#!/usr/bin/env python

import argparse
import numpy as np
from fasta import DNASequence

STOP, DIAG, UP, LEFT = 0, 1, 2, 3

def scoring_matrix(seq1: str, seq2: str, mode="global", match=1, mismatch=-1, gap=-1):
    """Claculates the scorint matrix for two candidate sequences for pairwise alignment
    :param seq1: str: string of nucleotide sequences for the first DNA sequence
    :param seq2: str: string of nucleotide sequences for the second DNA sequence
    :param mode: str: alignment node. options are local and global
    :param match int: score for a matching base between the two sequences, defalt is 1
    :param mismatch int: penalty for taking a mismatch to the next cell, default is -1
    :param gap int: penalty for taking a gap to the next cell, defaulu is -2
    """
    if mode not in {"global","local"}:
        raise ValueError("only valid entries for mode are global and local")
    score_mat = np.zeros((len(seq1)+1, len(seq2)+1), dtype=int)
    trace_mat = np.zeros((len(seq1)+1, len(seq2)+1), dtype=int)
    if mode == "global":
        for i in range(len(seq1)+1): #  O(n)
            score_mat[i,0]=i*gap     # O(1)
            trace_mat[i,0]= UP
        for j in range(len(seq2)+1): # O(m)
            score_mat[0,j]=j*gap     # O(1)
            trace_mat[0,j]=LEFT

        trace_mat[0, 0] = STOP
    for i in range(1, len(seq1)+1):   # O(n*m)
        for j in range(1, len(seq2)+1):
            diag = score_mat[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch) # O(1)
            up = score_mat[i-1, j] + gap   # O(1)
            left = score_mat[i, j-1] + gap  # O(1)


            if mode == "global":
                scores = [diag, up, left]
                best= max(scores)
                trace_mat[i, j] = [DIAG, UP, LEFT][scores.index(best)]
            else:
                scores = [0, diag, up, left]
                best= max(scores)
                trace_mat[i, j] = scores.index(best)
            score_mat[i, j] = best # O(1)

    return score_mat, trace_mat # mode

def back_track(seq1, seq2, score_mat, trace_mat, mode="global"):
    align1 = []
    align2 = []

    if mode == "global":
        i, j = len(seq1), len(seq2)
    elif mode == "local":
        i, j = np.unravel_index(np.argmax(score_mat), score_mat.shape)
    else:
        raise ValueError("mode must be 'global' or 'local'")


    while True:
        if mode == "global" and i == 0 and j == 0:
            break
        if mode == "local" and score_mat[i,j] == 0:
            break

        direction = trace_mat[i, j]

        if direction == DIAG:  # diagonal(match)
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1

        elif direction == UP:  #  up (gap in seq2)
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1

        elif direction == LEFT:  #  left (gap in seq1)
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1
        elif direction == STOP:
            break

    return "".join(reversed(align1)), "".join(reversed(align2))

def parse_args():
    parser = argparse.ArgumentParser(description="Sequence alignment tool for DNA sequences")
    parser.add_argument("ifile",help="Input file containing DNA sequnces in fasta format(.fa, .fas, .fasta)")
    subparsers = parser.add_subparsers(dest="mode",required=True)

    pairwise = subparsers.add_parser("pairwise",help="Pairwise Sequence alignment")
    pairwise.add_argument("--needle",action="store_true",help="Global alignment (Needleman-Wunsch)")
    pairwise.add_argument("--smith",action="store_true",help="Local alignment (Smith-Waterman)")
    pairwise.add_argument("--seq1",required=True,help="ID of first sequence")
    pairwise.add_argument("--seq2",required=True,help="ID of second sequence")

    msa = subparsers.add_parser("msa",help="Multiple sequence alignment")
    return parser.parse_args()

def main():
    print("Programm initiated")
    args = parse_args()
    #sequences = args.file
    sequences = DNASequence.fasta_file(args.ifile)

    if args.needle:
        seq1 = sequences[args.seq1].seq
        seq2 = sequences[args.seq2].seq
        s_mat, t_mat = scoring_matrix(seq1, seq2, mode="global")
        aln1, aln2 = back_track(seq1, seq2, s_mat, t_mat, mode="global")
    elif args.smith:
        seq1 = sequences[args.seq1].seq
        seq2 = sequences[args.seq2].seq
        s_mat, t_mat = scoring_matrix(seq1, seq2, mode="local")
        aln1, aln2 = back_track(seq1, seq2, s_mat, t_mat, mode="local")
    #print("Scoring Matrix:\n", s_mat)
    #print("\nTracking Matrix:\n", t_mat)

    print(aln1)
    print(aln2)
    print("Finished!")

if __name__ == "__main__":
    main()
