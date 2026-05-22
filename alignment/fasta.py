#!/usr/bin/env python3

import argparse
import re

class Sequence:
    """DNA sequence normalization and analysis utilities"""
    def __init__(self, seq):
        self.seq = self._normalize_sequence(seq)
        if not self.seq:
            raise ValueError("Empty DNA sequence")

    @property
    def gc_content(self):
        """Calculates the percentage of the DNA sequence compose og G and C nucleotides"""
        nuc_count = self.nucleotide_count()
        total = self.length
        if total == 0:
            raise ValueError("Cannot return gc content of empty sequence")
        gc_count = nuc_count['G'] + nuc_count['C']
        gc = gc_count/total
        return gc

    @property
    def length(self):
        return len(self.seq)

    @property
    def reverse_complement(self):
        complement = {'A':'T', 'T':'A','C':'G', 'G':'C'}
        return "".join(complement[base] for base in reversed(self.seq))

    @staticmethod
    def _normalize_sequence(seq):
        """Normalizes DNA/Protein sequence:
        - eliminates whitespaces and newline characters
        - uppercase all bases
        - validates nucleotides
        """
        valid_dna= {'A','C','G','T'}
        valid_prot = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'}
        seq=seq.upper()
        seq="".join(seq.split())   #seq.replace(" ", "").replace("\n", "")
        for i in seq:
            if i not in valid_dna and i not in valid_prot:
                raise ValueError("Invalid sequence!")
        return seq


    @staticmethod
    def fasta_parser(string):
        """
        parse FASTA formatted string
        :parameter string: string of DNA sequence in FASTA format
        :return sequences: dictionary of sequence IDs as keys and sequence string as values
        """

        sequences = {}
        current_id = None
        current_seq = []
        for line in string.splitlines():
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = Sequence("".join(current_seq))


                # Uniprot header
                m = re.match(r"^>sp\|[^|]+\|([^ ]+)", line)
                if m:
                    current_id = m.group(1)
                else:
                    current_id = line[1:] # fallback to full header without '>'

                current_seq  = []
            else:
                current_seq.append(line)

        if current_id is not None:
            sequences[current_id] = Sequence("".join(current_seq))

        return sequences

    @staticmethod
    def gc_contents(sequences):
        """Calculates GC content of multiple sequences
        :param *sequences: multiple ssequences in FASTA format
        :return gc_all: list of floats of gc content for all sequences
        """
        if not sequences:
            raise ValueError("No sequences provided")

        return [dna.gc_content for dna in sequences.values()]

    def nucleotide_count(self):
        """Counts number of nucleotides in a DNA sequences"""
        seq = self.seq
        if not seq:
            raise ValueError("Validation failed!")

        nuc_count = {'A':0,'C':0,'G':0,'T':0}
        for n in seq:
            nuc_count[n] += 1
        return nuc_count

    @staticmethod
    def fasta_file(path):
        if not path.lower().endswith((".fa",".fasta",".fna")):
            raise ValueError("file may not be valid")

        with open(path) as f:
            print("Read successful")
            return Sequence.fasta_parser(f.read())


def parse_args():
    parser = argparse.ArgumentParser(description="DNA/Protein sequence analysis from FASTA files")
    parser.add_argument("fasta_file", help="Input file, currently only FASTA format allowed (.fa, .fasta, .fna)")
    parser.add_argument("--output",help="Output file (default: stdout)", default=None)
    parser.add_argument("--gc",action="store_true",help="Report GC content")
    parser.add_argument("--rc",action="store_true",help="Report reverse compliment")

    return parser.parse_args()


def main():
    args = parse_args()
    sequences = Sequence.fasta_file(args.fasta_file)
    lines = []
    for seq_id, seq in sequences.items():
        #print(f"{seq_id}: {seq.seq} | Length: {seq.length}")
        fields = [seq_id, str(seq.length)]
        if args.gc:
            fields.append(f"{seq.gc_content:.2f}")
            #print(f"   GC: {seq.gc_content:.2f}")
        if args.rc:
            fields.append(f"{seq.reverse_complement}")
            #print(f"RC = {seq.reverse_complement}")

        lines.append("\t".join(fields))

    output_text = "\n".join(lines)

    if args.output:
        path = args.output
        with open(path, "w") as f:
            f.write(output_text + "\n")
            print("Write Successful")
    else:
        print(output_text)

    #print(gc_content(seq))
    # argument used in function call must match function definition. omitting the star results in a nested tuple



if __name__ == "__main__":
    main()
