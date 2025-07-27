# Needleman-Wunsch Pairwise Alignment Tool

This project implements a command-line application to perform global pairwise sequence alignment using the Needleman-Wunsch algorithm. It takes two input sequences from a `.fna` file and a scoring matrix from a `.mtx` file, performs alignment, and outputs the result in FASTA-like format.

---

## Features

- Parses `.fna` files containing biological sequences in FASTA format
- Parses `.mtx` files containing scoring matrices (e.g., nucleotide or BLOSUM)
- Implements the Needleman-Wunsch algorithm for global sequence alignment
- Outputs aligned sequences and score to a user-defined `.fna` file

---

## File Structure

| File | Description |
|------|-------------|
| `nw_solver.py`     |  **Main entry point** — orchestrates I/O and alignment |
| `cli_auth.py`      | Validates and parses command-line arguments |
| `fna_tools.py`     | Reads `.fna` files and writes aligned output |
| `mtx_tools.py`     | Parses scoring matrices and handles score lookup |
| `nw_alignment.py`  | Performs matrix initialization, alignment scoring, and traceback |

---

## Dependencies

- Python 3.7+
- `numpy`

Install with:

```bash
pip install numpy
