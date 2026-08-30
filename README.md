# Sequence Aligner Project

Implemented pairwise sequence alignment using dynamic programming, with fasta file parsing. 

This project provides a command-line interface tool to parse FASTA files, analyze sequence statistics (length, GC content, and reverse complements), and run alignments using either standard simple gap penalties or affine gap penalties (opening/extension) with optional PAM and BLOSUM substitution matrices.

> [!NOTE]
> *Pending: Alignment with HMM (Hidden Markov Model).*

---

## Project Structure

*   **`alignment/`** — Core alignment scripts and data files.
    *   [`fasta.py`](file:///c:/Users/sidne/OneDrive/Desktop/personal/Sequence-Aligner_project/alignment/fasta.py) — Parser for FASTA formatted sequence files and analysis utility.
    *   [`aln.py`](file:///c:/Users/sidne/OneDrive/Desktop/personal/Sequence-Aligner_project/alignment/aln.py) — Dynamic programming pairwise alignment (Needleman-Wunsch / Smith-Waterman).
    *   [`affine_align.py`](file:///c:/Users/sidne/OneDrive/Desktop/personal/Sequence-Aligner_project/alignment/affine_align.py) — Affine gap penalty alignment using PAM/BLOSUM.
    *   [`PAM250.txt`](file:///c:/Users/sidne/OneDrive/Desktop/personal/Sequence-Aligner_project/alignment/PAM250.txt) & [`Blosum62.txt`](file:///c:/Users/sidne/OneDrive/Desktop/personal/Sequence-Aligner_project/alignment/Blosum62.txt) — Substitution matrices.
*   **`FasraParser/`** — Playground directory containing test sequences and early prototypes (`w1.py`, `w3.py`).

---

## Requirements

The project relies on standard Python libraries and **NumPy** for scoring/traceback matrices.
```bash
pip install numpy
```

---

## Usage & Examples

### 1. Analyzing FASTA Files (`fasta.py`)
Extract basic properties such as sequence lengths, GC contents, and reverse complements.

```bash
python alignment/fasta.py alignment/sample_a.fasta --gc --rc
```

**Output:**
```text
Read successful
seqA_001	6	0.50	ACGCAT
seqA_002	7	0.29	TGTAATC
seqA_003	6	0.50	TAACGG
...
seqA_012	6	0.00	ATATAT
```

---

### 2. Pairwise Dynamic Programming Alignment (`aln.py`)
Perform global (Needleman-Wunsch) or local (Smith-Waterman) alignments.

#### Global Alignment Example:
```bash
python alignment/aln.py alignment/sample_a.fasta pairwise --needle --seq1 seqA_001 --seq2 seqA_002
```
**Output:**
```text
Programm initiated
Read successful
-A-TGCGT
GATTAC-A
Finished!
```

#### Local Alignment Example:
```bash
python alignment/aln.py alignment/sample_a.fasta pairwise --smith --seq1 seqA_001 --seq2 seqA_002
```
**Output:**
```text
Programm initiated
Read successful
AT
AT
Finished!
```

---

### 3. Alignment with Affine Gap Penalties (`affine_align.py`)
Align sequences using separate gap opening (`go`) and gap extension (`ge`) penalties. You can also specify PAM or BLOSUM matrices.

#### Affine Alignment Example:
```bash
python alignment/affine_align.py alignment/sample_a.fasta pairwise --needle --seq1 seqA_001 --seq2 seqA_002
```
**Output:**
```text
Read successful
Initiated
-ATGCGT
GATTACA
Finished!
```

#### Affine Alignment with PAM250 Matrix Example:
```bash
python alignment/affine_align.py alignment/sample_a.fasta pairwise --needle --seq1 seqA_001 --seq2 seqA_002 --dna alignment/PAM250.txt
```
**Output:**
```text
Read successful
Initiated
--ATGCGT
GATTAC-A
Finished!
```

---

## Algorithmic Details

*   **Time Complexity**:
    *   Scoring matrix calculation runs in $\mathcal{O}(n \times m)$ time, where $n$ and $m$ are the lengths of the two sequences.
    *   Sequence normalization and ID lookup run in linear time $\mathcal{O}(n)$ relative to sequence length.
*   **Traceback States**:
    *   For standard alignment: `DIAG` (diagonal match/mismatch), `UP` (gap in sequence 2), and `LEFT` (gap in sequence 1).
    *   For affine alignment: Three matrices ($M$, $X$, $Y$) track matches, vertical gap transitions, and horizontal gap transitions respectively.
