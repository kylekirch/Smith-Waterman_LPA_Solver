# Smith-Waterman Local Pairwise Alignment (LPA)

This repository implements a command-line tool for **local sequence alignment** using the **Smith-Waterman algorithm**. It takes biological sequences in `.fna` format and a scoring matrix in `.mtx` format, computes the optimal local alignment, and outputs the aligned sequences with their alignment score in a `.fna` file.

---

## Features

- **Local alignment** via the Smith-Waterman algorithm
- Parses input sequences from FASTA-like `.fna` files
- Parses scoring matrices from `.mtx` files
- Produces output in FASTA-like `.fna` format, with alignment scores
- CLI argument validation for input, output, and scoring files
- Modular design (`cli_auth`, `fna_tools`, `mtx_tools`, `sw_alignment`, `sw_solver`)

---

## File Overview

| File            | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `sw_solver.py`  | **Main entry point** — orchestrates input/output and alignment process       |
| `cli_auth.py`   | Handles and validates CLI arguments                                         |
| `fna_tools.py`  | Parses `.fna` input sequences and writes aligned sequences to output         |
| `mtx_tools.py`  | Parses `.mtx` scoring matrix files and provides lookup utilities             |
| `sw_alignment.py` | Contains Smith-Waterman logic: matrix initialization, scoring, traceback |

---

## Requirements

- Python 3.7+
- NumPy

Install dependencies with:

```bash
pip install numpy
