"""
Filename: sw_alignment.py
Project: Smith-Waterman Local Pairwise Alignment (LPA)
Description: Implements the Smith-Waterman algorithm, including
              score matrix initialization, scoring, and traceback logic
              for local pairwise sequence alignment.
Author: Kyle Kirchgessner
Date: 2025-07-28
Version: 1.0
Dependencies: numpy, mtx_tools
"""


from mtx_tools import query_mtx
import numpy as np


def sw_initialize_matrices(num_rows, num_cols):
    '''
    Description:
    Initializes the Smith-Waterman score and traceback matrices.

    IN:
    para1: desired number of rows; [int]
    para2: desired number of columns; [int]
    OUT:
    return1: score matrix filled with 0s; [np.ndarray]
    return2: traceback matrix filled with 'X'; [2D list]
    '''

    # fill first row with default values
    score_matrix = [[0] * (num_cols + 1)]
    traceback_matrix = [['X'] * (num_cols + 1)]

    # fill first column with default values
    for _ in range(num_rows):
        score_matrix.append([0])
        traceback_matrix.append(['X'])

    return score_matrix, traceback_matrix


def sw_traceback(sequence_1, sequence_2, top_score_location, traceback_matrix):
    '''
    Description:
    Performs a needleman-wunsch traceback

    IN:
    para1: first sequence; [string]
    para2: second sequence; [string]
    para3: traceback matrix; [2d npArray or list of lists]

    OUT:
    return1: the traced sequences (sequence_1_aligned, sequence_2_aligned); [tuple(2)]
    '''

    # initialize the starting position as the position of the largest score
    row = top_score_location[0]
    col = top_score_location[1]

    # initialize the final paired sequences
    sequence_1_aligned = ""
    sequence_2_aligned = ""

    # perform traceback steps until [0][0] is reached
    while True:
        # get the current direction from the traceback matrix
        direction = traceback_matrix[row][col]

        if direction == 'D':
            # match/mismatch: add both characters
            sequence_1_aligned += sequence_1[col - 1]
            sequence_2_aligned += sequence_2[row - 1]
            row -= 1
            col -= 1  # move up and left
        elif direction == 'U':
            # gap in sequence_1
            sequence_1_aligned += '-'
            sequence_2_aligned += sequence_2[row - 1]
            row -= 1  # move up
        elif direction == 'L':
            # gap in sequence_2
            sequence_1_aligned += sequence_1[col - 1]
            sequence_2_aligned += '-'
            col -= 1  # move left
        elif direction == 'X':
            break
        else:
            # if anything other than the allowed traces are detected, throw error
            raise TypeError("NW_TRACEBACK ERROR: UNKNOWN TRACE DETECTED")

        # exit once [0][0] is reached
        if row == 0 and col == 0:
            break

    # return the paired sequences in the correct orientation
    return (sequence_1_aligned[::-1], sequence_2_aligned[::-1])


def sw_alignment(sequence_1, sequence_2, score_matrix):
    '''
    Description:
    Performs needleman-wunsch pairwise alignment

    IN:
    para1: the first sequence; [string]
    para2: the second sequence; [string]
    para3: the scoring matrix; [2d npArray]

    OUT:
    return1: the alignment score; [int]
    return2: the traced sequences (sequence_1_aligned,sequence_2_aligned): [tuple(2)]
    '''

# ensure the existance of a global score table (eg. BLOSUM50, nucleotide)
    if (score_matrix is not None):
        score_table = score_matrix
        global SCORING_MATRIX
        SCORING_MATRIX = score_table
    else:  # if no score table exists, throw error
        raise ValueError('NW_ALIGNMENT ERROR: No score matrix provided')

  # initial parameters
    num_cols = len(sequence_1)
    num_rows = len(sequence_2)
    gap = query_mtx('A', '-')
    final_score_matrix, traceback_matrix = sw_initialize_matrices(
        num_rows, num_cols)

  # initialize a top score
    top_score = [0, (0, 0)]

  # perform needleman wunsch pairwise alignment
    for row in range(num_rows):
        row = row+1  # ensure the default value is skipped
        for col in range(num_cols):
            col = col+1
            possible_scores = {
                'D': final_score_matrix[row-1][col-1] + query_mtx(sequence_1[col-1], sequence_2[row-1]),
                'L': final_score_matrix[row][col-1] + gap,
                'U': final_score_matrix[row-1][col] + gap,
            }
            # determine which option has the maximum value, store the value in the scores table, and store the trace(s) in the traceback table
            max_val = max(possible_scores.values())

            # if the max value is less than 0, use 0 as the score, and 'X' as the trace
            if (max_val <= 0):
                max_val = 0
                trace = 'X'
                # add the max score and trace to their respective matrices
                final_score_matrix[row].append(max_val)
                traceback_matrix[row].append(trace)
                continue
            if (max_val >= top_score[0]):
                top_score = (max_val, (row, col))
            for name, val in possible_scores.items():
                if val == max_val:
                    trace = name  # in the case of multiple matching max values, it prioritizes the lowest index in the hashtable, in this case 'D'
                    break

            # add the max score and trace to their respective matrices
            final_score_matrix[row].append(max_val)
            traceback_matrix[row].append(trace)

    scoreArray = np.array(final_score_matrix)
    tracebackArray = np.array(traceback_matrix)
    top_score_location = top_score[1]

    sequence_1_final, sequence_2_final = sw_traceback(
        sequence_1, sequence_2, top_score_location, tracebackArray)
    return (top_score[0], sequence_1_final, sequence_2_final)


# debug
ROW_LABELS = ['A', 'C', 'G', 'T', '-']
COL_LABELS = ['A', 'C', 'G', 'T', '-']
SCORING_MATRIX = [[2, -1, -1, -1, -2],
                  [-1, 2, -1, -1, -2],
                  [-1, -1, 2, -1, -2],
                  [-1, -1, -1, 2, -2],
                  [-2, -2, -2, -2, 2]]


def main():

    sequence_1 = 'AGC'
    sequence_2 = 'AAAC'
    top_score, top_score_location, sequence_1_final, sequence_2_final = sw_alignment(
        sequence_1, sequence_2, SCORING_MATRIX)
    print(top_score, top_score_location, sequence_1_final, sequence_2_final)


if __name__ == "__main__":
    main()
