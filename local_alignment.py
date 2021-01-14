import numpy as np

# Local alignment of 2 sequences using the Smith-Waterman Algorithm



def compute_local_alignment(sequence1, sequence2, gap , match, mismatch):
    l1 = len(sequence1)
    l2 = len(sequence2)


    print('Local Alignment')

    # Step 1: Instantiate the matrix
    smith_wat_matrix = instantiate_matrix(l1, l2)

    # Step 2: Compute matrix values
    smith_wat_matrix = compute_matrix(smith_wat_matrix, sequence1, sequence2, gap , match , mismatch)

    print('Scoring Matrix:')
    print(smith_wat_matrix)
    # Step 3: Compute boolean traceback matrix
    # Finds the indexes of the maximum element of the matrix
    ind = np.unravel_index(np.argmax(smith_wat_matrix, axis=None), smith_wat_matrix.shape)
    trace_back = compute_traceback_matrix(smith_wat_matrix, 0, 0, None, ind)
    print('Alignment Path:')
    print(trace_back)

    # Step 4: Compute final alignment
    alignment1, alignment2 = compute_alignment(trace_back, sequence1, sequence2)

    print('Alignment For Sequences')
    print('(', sequence1, ',', sequence2, ')')
    alignment = (alignment1, alignment2)
    print('is:')
    return alignment


def instantiate_matrix(l1, l2):
    smith_wat_matrix = np.zeros((l1 + 1, l2 + 1), dtype='float')

    for index in range(1, l2 + 1):
        smith_wat_matrix[0, index] = 0

    for index in range(1, l1 + 1):
        smith_wat_matrix[index, 0] = 0
    return smith_wat_matrix


def compute_matrix(matrix, sequence1, sequence2, gap_penalty, match_penalty, mismatch_penalty):
    # Fills the score matrix
    for i in range(1, np.shape(matrix)[0]):
        for j in range(1, np.shape(matrix)[1]):
            row = matrix[i][j - 1] - gap_penalty
            col = matrix[i - 1][j] - gap_penalty
            # match
            if sequence1[i - 1] == sequence2[j - 1]:
                diag = matrix[i - 1][j - 1] + match_penalty
            # mismatch
            else:
                diag = matrix[i - 1][j - 1] - mismatch_penalty

            if max(row, col, diag) <= 0:
                matrix[i][j] = 0
            else:
                matrix[i][j] = (max(row, col, diag))

    return matrix


def compute_traceback_matrix(matrix, col_index, row_index, trace, max_ind):
    # Compute index of max element
    l1 = (np.shape(matrix))[0]
    l2 = (np.shape(matrix))[1]

    if col_index == 0 and row_index == 0:
        trace = np.zeros((l1, l2), dtype=bool)
        trace[max_ind[0], max_ind[1]] = True

        # Compute element on top , left and diagonal
        col = matrix[max_ind[0] - 1][max_ind[1]]
        row = matrix[max_ind[0]][max_ind[1] - 1]
        diag = matrix[max_ind[0] - 1][max_ind[1] - 1]
        # Finds the max and set it to true
        if diag >= col and diag >= row:
            trace[max_ind[0] - 1][max_ind[1] - 1] = True
            compute_traceback_matrix(matrix, col_index + 1, row_index + 1, trace, max_ind)
        elif col > diag and col >= row:
            trace[max_ind[0] - 1][max_ind[1]] = True
            compute_traceback_matrix(matrix, col_index + 1, row_index, trace, max_ind)
        else:
            trace[max_ind[0]][max_ind[1] - 1] = True
            compute_traceback_matrix(matrix, col_index, row_index + 1, trace, max_ind)
    # Recursive call
    else:
        if col_index < max_ind[0] and row_index < max_ind[1]:
            col = matrix[max_ind[0] - col_index - 1][max_ind[1] - row_index]
            row = matrix[max_ind[0] - col_index][max_ind[1] - 1 - row_index]
            diag = matrix[max_ind[0] - col_index - 1][max_ind[1] - 1 - row_index]

            if col == 0 or row == 0 or diag == 0:
                return

            # Finds the max and set it to true
            if diag >= col and diag >= row:
                trace[max_ind[0] - col_index - 1][max_ind[1] - 1 - row_index] = True
                compute_traceback_matrix(matrix, col_index + 1, row_index + 1, trace, max_ind)
            elif col > diag and col >= row:
                trace[max_ind[0] - col_index - 1][max_ind[1] - row_index] = True
                compute_traceback_matrix(matrix, col_index + 1, row_index, trace, max_ind)
            else:
                trace[max_ind[0] - col_index][max_ind[1] - 1 - row_index] = True
                compute_traceback_matrix(matrix, col_index, row_index + 1, trace, max_ind)
    return trace


def compute_alignment(trace, sequence1, sequence2):
    seq1 = ''
    seq2 = ''
    hist_i = 0
    hist_j = 0
    start_index = True

    for i in range(len(sequence1)):
        for j in range(len(sequence2)):
            if trace[i][j]:

                if start_index:
                    print(seq1)
                    seq1 += sequence1[i - 1]
                    seq2 += sequence2[j - 1]
                    start_index = False
                elif i == hist_i + 1 and j == hist_j + 1:

                    seq1 += sequence1[i - 1]
                    seq2 += sequence2[j - 1]

                elif i == hist_i and j == hist_j + 1:
                    seq1 += '-'
                    seq2 += sequence2[hist_i]
                elif i == hist_i + 1 and j == hist_j:
                    seq1 += sequence1[hist_j]
                    seq2 += '-'
                hist_i = i
                hist_j = j

    return seq1, seq2
