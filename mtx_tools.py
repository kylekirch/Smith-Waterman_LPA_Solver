import numpy as np
ROW_LABELS = []
COL_LABELS = []
SCORING_MATRIX = [[]]


def parse_mtx(filepath):
    '''
    Description: 
    Parses a given .mtx file and returns its indexable contents

    IN:
    para1: path to .mtx file; [string]
    .mtx file can contain an matrix of any rectangular dimension

    OUT:
    return1: inner matrix data; [2d npArray]
    return2: row labels; [list]
    return3: column labels; [list]
    '''

    with open(filepath, 'r') as file:
        lines = file.readlines()

    # store the first row as column labels
    col_labels = lines[0].strip().split()

    row_labels = []
    data = []

    # parse file line by line, skipping the row label and storing it
    for line in lines[1:]:
        entries = line.strip().split()
        if not entries:
            continue
        row_labels.append(entries[0])
        data.append([int(n) for n in entries[1:]])

    # convert to npArray, typecast data to integer type
    matrix = np.array(data, dtype=int)

    global ROW_LABELS
    global COL_LABELS
    global SCORING_MATRIX
    ROW_LABELS = row_labels
    COL_LABELS = col_labels
    SCORING_MATRIX = matrix

    return matrix, row_labels, col_labels

# debug (simulate global scoring matrix and row/column labels)


def query_mtx(target_row, target_column, matrix=None, rows=None, cols=None):
    '''
    Description:
    Returns a value from the provided matrix at the given position (target_row, target_column)

    IN:
    para1: the label of the target row; [string]
    para2: the label of the target column; [string]
    para3: an optional non-global scoring matrix; [2d npArray]
    para4: an optional non-global list of row labels [list];
    para5: an optional non-global list of column labels [list];

    OUT:
    return1: the value at the target position; [int]
    '''

    # if using a non-global scoring matrix
    if (matrix is not None):
        if (rows is not None and cols is not None):
            i = rows.index(target_row)
            j = cols.index(target_column)
            return matrix[i][j]
        else:
            raise ValueError(
                'ERROR: Custom matrix was provided without matching column/row labels')

    # index the previously parsed row and column labels
    i = ROW_LABELS.index(target_row)
    j = COL_LABELS.index(target_column)

    # return score from global matrix
    return SCORING_MATRIX[i][j]


def main():
    row_labels = ['A', 'C', 'G', 'T', '-']
    col_labels = ['-', 'C', 'T', 'A', 'G']
    matx = [[3, -3, -3, -3, -6],
            [-3, 3, -3, -3, -6],
            [-3, -3, 3, -3, -6],
            [-3, -3, -3, 3, -6]]

    result1 = query_mtx('A', 'T')
    result2 = query_mtx('C', 'C')
    result3 = query_mtx('-', 'G')

    print(result1, result2, result3, '\n')

    result1 = query_mtx('A', 'T', matx, row_labels, col_labels)
    result2 = query_mtx('C', 'C', matx, row_labels, col_labels)
    result3 = query_mtx('G', 'T', matx, row_labels, col_labels)

    print(result1, result2, result3, '\n')


if __name__ == "__main__":
    main()
