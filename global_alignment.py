import numpy as np

# Global alignment of 2 sequences using the Needleman-Wunsch Algorithm

gap_penalty = 1
match_penalty = 1


def compute_global_alignment(sequence1, sequence2, gap , match , mismatch ):
    alignment = ''
    l1 = len(sequence1)
    l2 = len(sequence2)

    print('Global Alignment')
    # Step 1: Instantiate the matrix
    needleMatrix = instantiate_matrix(l1, l2)

    # Step 2: Fill in the table
    needleMatrix = compute_matrix(needleMatrix, sequence1, sequence2, gap, match , mismatch)
    print('Scoring Matrix')
    print(needleMatrix)
    # Step 3: Trace-Back
    trace_back_matrix = trace_back(needleMatrix, 0, 0, None)
    print('Path Matrix')
    print(trace_back_matrix)
    # Step 4: Process Trace-Back
    alignment1, alignment2 = compute_alignment(trace_back_matrix, sequence1, sequence2)
    print('Alignment For Sequences')
    print('(', sequence1, ',', sequence2, ')')
    print('Is:')
    alignment = (alignment1, alignment2)
    return alignment


def instantiate_matrix(l1, l2):
    needleMatrix = np.zeros((l1 + 1, l2 + 1), dtype='float')

    for index in range(1, l2 + 1):
        needleMatrix[0, index] = -index

    for index in range(1, l1 + 1):
        needleMatrix[index, 0] = -index
    return needleMatrix


def compute_matrix(matrix, sequence1, sequence2, gap_penalty, match_penalty, mismatch_penalty):
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
            matrix[i][j] = (max(row, col, diag))

    return matrix


# Recursive function
def trace_back(matrix, col_index, row_index, trace):
    l1 = (np.shape(matrix))[0]
    l2 = (np.shape(matrix))[1]
    # Reaches the end of the matrix
    if col_index + 1 == l1 or row_index + 1 == l2:
        return trace

    # Base case
    if col_index == 0 and row_index == 0:
        trace = np.zeros((l1, l2), dtype=bool)
        trace[l1 - 1, l2 - 1] = True

        col = matrix[l1 - 2][l2 - 1]
        row = matrix[l1 - 1][l2 - 2]
        diag = matrix[l1 - 2][l2 - 2]

        if diag >= col and diag >= row:
            trace[l1 - 2][l2 - 2] = True
            trace_back(matrix, col_index + 1, row_index + 1, trace)
        elif col > diag and col >= row:
            trace[l1 - 2][l2 - 1] = True
            trace_back(matrix, col_index + 1, row_index, trace)
        else:
            trace[l1 - 1][l2 - 2] = True
            trace_back(matrix, col_index, row_index + 1, trace)

    # Recursive call
    else:
        col = matrix[l1 - col_index - 1][l2 - row_index]
        row = matrix[l1 - col_index][l2 - row_index - 1]
        diag = matrix[l1 - col_index - 1][l2 - row_index - 1]

        if diag >= col and diag >= row:
            trace[l1 - col_index - 1][l2 - row_index - 1] = True
            trace_back(matrix, col_index + 1, row_index + 1, trace)
        elif col > diag and col >= row:
            trace[l1 - col_index - 1][l2 - row_index] = True
            trace_back(matrix, col_index + 1, row_index, trace)
        else:
            trace[l1 - col_index][l2 - row_index - 1] = True
            trace_back(matrix, col_index, row_index + 1, trace)

    return trace


# Given the boolean matrix , try to find the best alignment possible
def compute_alignment(trace, sequence1, sequence2):
    seq1 = ''
    seq2 = ''
    hist_i = 0
    hist_j = 0
    start_index = True

    for i in range(len(sequence1)+1):
        for j in range(len(sequence2)+1):
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
                    seq2 += sequence2[j-1]
                elif i == hist_i + 1 and j == hist_j:
                    seq1 += sequence1[i-1]
                    seq2 += '-'
                hist_i = i
                hist_j = j

    return seq1, seq2
