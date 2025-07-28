"""
Filename: sw_solver.py
Project: Smith-Waterman Local Pairwise Alignment (LPA)
Description: Main entry point for the application. Orchestrates input parsing,
              sequence alignment using Smith-Waterman, and output generation.
Author: Kyle Kirchgessner
Date: 2025-07-28
Version: 1.0
Dependencies: cli_auth, fna_tools, mtx_tools, sw_alignment
"""

from cli_auth import cli_auth
import fna_tools as fna
import mtx_tools as mtx
import sw_alignment as sw


def main():
  # declare global variables
    global SCORING_MATRIX
    global COL_LABELS
    global ROW_LABELS

  # extract path arguments from command line
    args = cli_auth()
    input_file = args.i
    output_file = args.o
    score_file = args.s

  # extract sequence and scoring matrix data from their respective files
    sequence_data = fna.parse_fna(input_file)
    SCORING_MATRIX, ROW_LABELS, COL_LABELS = mtx.parse_mtx(score_file)
    sequence_1 = sequence_data[0][2]
    sequence_2 = sequence_data[1][2]

  # perform smith-waterman local alignment
    top_score, sequence_1_aligned, sequence_2_aligned = sw.sw_alignment(
        sequence_1, sequence_2, SCORING_MATRIX)
    print('\n', sequence_1_aligned, '\n', sequence_2_aligned, '\n', top_score)

  # print results to output file
    sequence_1_id = sequence_data[0][0]
    sequence_2_id = sequence_data[1][0]
    fna.write_pairs(output_file, (sequence_1_id, sequence_1_aligned),
                    (sequence_2_id, sequence_2_aligned), top_score)


if __name__ == "__main__":
    main()
