n = 22
r = 14
k = n - r
l = 3
Fq = GF(2)

"""Perform Partial Gaussian Elimination (PGE) on given matrix.

Args:
    A: matrix to perform PGE on
    l: PGE argument

Returns:
    If PGE can be applied, returns change of basis matrix S.
    Otherwise returns False.
"""
def pge(A, l):

    #Check rank to see if it's possible to apply PGE
    total_rank = A[:, -(r-l):].rank()
    if total_rank < r-l:
        return False

    # Change of basis matrix initialization
    S = identity_matrix(Fq, r)
    # Column of the element to pivot on
    pivot_column = k+l

    # Pivot on each row
    for pivot_row in range(r-l):

        pivot = A[pivot_row][pivot_column]
        row_to_check = pivot_row+1

        # If the pivot is null, find a row to swap
        while pivot == 0 and row_to_check < r:
            if A[row_to_check][pivot_column] == 1:
                P = identity_matrix(Fq, r) 
                P.swap_rows(pivot_row, row_to_check)
                A = P*A
                S = P*S
                pivot=A[pivot_row][pivot_column]
            row_to_check += 1

        # Null all the other elements on the same column
        for j in range(r):
            if j!=pivot_row and A[j][pivot_column] == 1:
                A = sum_rows(A, j, pivot_row)
                M = identity_matrix(Fq, r) 
                M[j, pivot_row] = 1
                S = M*S

        # Increment column index
        pivot_column+=1


    # If it's partial GE, rows need to be swapped to obtain desired form
    if l != 0:
        P = identity_matrix(Fq, r)
        P = P[-l:,:].stack(P[:r-l,:])
        A = P*A
        S = P*S

    # Return the change of basis matrix
    return S


"""Perform sum of rows on a matrix

Args:
    A: matrix to perform PGE on
    affected_row: row that will change
    addend_row: row that is summed to the affected row

Returns:
    Return the changed matrix
"""
def sum_rows(A, affected_row, addend_row):
    for i in range(A.ncols()):
        A[affected_row, i] = A[affected_row][i] + A[addend_row][i]
    return A
        
A = random_matrix(Fq, r, n)
B = matrix(A)
S = pge(A, l) #S <-- change of basis matrix
if S is not False:
    print(S*B)
else:
    print(f"Calculated A not good, rank of last r-l columns is {A[:, -(r-l):].rank()}")
