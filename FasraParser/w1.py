#!/usr/bin/env python3



class DNASequence:
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
        """Normalizes DNA sequence:
        - eliminates whitespaces and newline characters
        - uppercase all bases
        - validates nucleotides
        """
        valid= {'A','C','G','T'}
        seq=seq.upper()
        seq="".join(seq.split())   #seq.replace(" ", "").replace("\n", "")
        for i in seq:
            if i not in valid:
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
                    sequences[current_id] = DNASequence("".join(current_seq))
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            sequences[current_id] = DNASequence("".join(current_seq))

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
            return DNASequence.fasta_parser(f.read())

    def file_output(self, path='output.txt'):
        """Allows saving results to an output file of specified format"""
        pass

    def dummy_func(self, path):
        print(f"running {path}")

def usage():
    """Usage instructions for the python script
    Will adjust later to use command line arguments
    """
    print("This script will be used like this")


def main():
    #seq = input("Enter seq: ")
    fasta = """>seq1
    ATGC
    GTAA

    >seq2
    GGG
    """
    sequences= DNASequence.fasta_parser(fasta)
    for key, value in sequences.items():
        print(f"{key}: {value.seq} | GC: {value.gc_content:.2f}")
        print(f"RC = {value.reverse_complement}")
    dna = DNASequence("ATGC")
    dna.dummy_func("test.fasta")
    #print(gc_content(seq))
    # argument used in function call must match function definition. omitting the star results in a nested tuple



if __name__ == "__main__":
    main()
