from mtx_tools import query_mtx
import numpy as np


def nw_initialize_matrices(num_rows, num_cols, gap):
    '''
    Description:
    Initializes the Needleman-Wunsch final score matrix and traceback matrix

    IN:
    para1: desired number of rows; [int]
    para2: desired number of columns; [int]
    para3: gap penalty [int]

    OUT:
    return1: empty final score matrix, with only the first column and row filled with the default values [2d list]
    return2: empty traaceback, with only the first column and row filled with the default traces [2d list]
    '''

  # initialize the final scoring matrix
    max_dim = max(num_cols, num_rows)
    default_vals = [0]
    default_score_matrix = []
    default_traceback_matrix = []
  # determine the set of possible default values that will fill the first row and first column
    for index in range(max_dim):
        default_vals.append(default_vals[index] + gap)

  # initialize the first row of the scoring/traceback matrix with default values
    default_score_matrix = [default_vals[:num_cols+1]]
    default_traceback_matrix = [['D'] + ['L'] * (num_cols)]

  # initialize the first element of each row with its default value
    for i in range(num_rows):
        default_score_matrix.append([default_vals[i+1]])
        default_traceback_matrix.append(['U'])

    return (default_score_matrix, default_traceback_matrix)


def nw_traceback(sequence_1, sequence_2, traceback_matrix):
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

    # initialize the starting position as the last element in the traceback matrix
    num_cols = len(sequence_1)
    num_rows = len(sequence_2)
    row = num_rows
    col = num_cols

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
        else:
            # if anything other than the allowed traces are detected, throw error
            raise TypeError("NW_TRACEBACK ERROR: UNKNOWN TRACE DETECTED")

        # exit once [0][0] is reached
        if row == 0 and col == 0:
            break

    # return the paired sequences in the correct orientation
    return (sequence_1_aligned[::-1], sequence_2_aligned[::-1])


def nw_alignment(sequence_1, sequence_2, score_matrix):
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
        global SCORE_MATRIX
        SCORE_MATRIX = score_table
    else:  # if no score table exists, throw error
        raise ValueError('NW_ALIGNMENT ERROR: No score matrix provided')

  # initial parameters
    num_cols = len(sequence_1)
    num_rows = len(sequence_2)
    gap = query_mtx('A', '-')
    final_score_matrix, traceback_matrix = nw_initialize_matrices(
        num_rows, num_cols, gap)

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
            for name, val in possible_scores.items():
                if val == max_val:
                    trace = name  # in the case of multiple matching max values, it prioritizes the lowest index in the hashtable, in this case 'D'
                    break

            # add the max score and trace to their respective matrices
            final_score_matrix[row].append(max_val)
            traceback_matrix[row].append(trace)

    scoreArray = np.array(final_score_matrix)
    tracebackArray = np.array(traceback_matrix)
    sequence_1_final, sequence_2_final = nw_traceback(
        sequence_1, sequence_2, tracebackArray)
    score = scoreArray[-1, -1]
    return (score, sequence_1_final, sequence_2_final)


def main():
    score_matrix = [[1, -1, -1, -1, -2],
                    [-1, 1, -1, -1, -2],
                    [-1, -1, 1, -1, -2],
                    [-1, -1, -1, 1, -2],
                    [-2, -2, -2, -2, 1]]
    sequence_1 = 'AGC'
    sequence_2 = 'AAAC'

    score, sequence_1_paired, sequence_2_paired = nw_alignment(
        sequence_1, sequence_2, score_matrix)
    print('\n', sequence_1_paired, '\n', sequence_2_paired, '\n', score)


if __name__ == "__main__":
    main()
