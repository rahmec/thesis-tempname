import random
from sage.coding.code_constructions import random_linear_code

# Remember to adjust notation to the inernational pyhon standard before goin onward

# Define the various sizes

p0 = 1
n = 22
r = 11
k = n - r
l = 11
Fq = GF(2)

def generate_random_linear_code(domain, n, k):
    C = random_linear_code(domain, n, k)
    return C.parity_check_matrix()

def generate_random_ldpc(n,r,wr,wc):
    return generate_matrix(r,n,wc,wr)

'''
Args:
    m : number of rows
    n : number of columns
    wr: ones per row
    wc: ones per columns

It is mandatory that m*wr = n*wc (
'''

def generate_matrix(m, n, wr, wc):
    H = Matrix(m,n)
    total_ones_indices = [[] for _ in range(m)]
    for i in range(n):
        available_indices = [x for x in range(m)]
        for _ in range(wc):
            row = random.choice(available_indices)
            H[row, i]=1
            available_indices.remove(row)
            total_ones_indices[row].append(i)
    
    underweighted_rows = []
    overweighted_rows = []
    for i in range(m):
        #if len(total_ones_indices[i]) < 2:
        #    underweighted_rows.append(i)
        #    zero_indices = [r for r in range(n) if r not in total_ones_indices[i]]
        #    for _ in range(2-len(total_ones_indices[i])):
        #        col = random.choice(zero_indices)
        #        zero_indices.remove(col)
        #        H[i,col]=1
        #        total_ones_indices[i].append(col)
        if len(total_ones_indices[i]) < wr:
            underweighted_rows.append(i)
        elif len(total_ones_indices[i]) > wr:
            overweighted_rows.append(i)

    print(overweighted_rows)
    print(underweighted_rows)

    for ov_row in overweighted_rows:
        for _ in range(len(total_ones_indices[ov_row])-wr):
            found_flag = False
            available_indices = total_ones_indices[ov_row]
            while not found_flag:
                col = random.choice(available_indices)
                #print(H)
                for un_row in underweighted_rows:
                    if col not in total_ones_indices[un_row]:
                        H[un_row, col] = 1
                        H[ov_row, col] = 0
                        total_ones_indices[un_row].append(col)
                        total_ones_indices[ov_row].remove(col)
                        if len(total_ones_indices[un_row]) == wr:
                            underweighted_rows.remove(un_row)
                        found_flag = True
                        break
                if col in available_indices:
                    available_indices.remove(col)
    
    print(H)
    return H

'''Converts the permutation vector sigma into a permutation matrix
E.g. 
                1 0 0 0
                0 0 1 0 
[1 3 2 4] --->  0 1 0 0
                0 0 0 1

'''
def permutation_matrix(sigma):
    n = len(sigma)
    P = matrix(Fq, n)
    for i in range(n):
        P[i, sigma[i] - 1] = 1
    return P

'''Generates a random permutation vector of size n.
E.g. [1 3 2 4]
'''
def sample_permutation_matrix(n):
    sigma = Permutations(n).random_element()
    return permutation_matrix(sigma)

'''
ahecks if the leftmost k+l columns form a submatrix with rank >= k
'''
def check_random_permutation(PH, l):
    submatrix = PH[:, -(r-l):]
    return submatrix.rank() >= r-l

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

    # Return the changed matrix
    return A

def lee_brickell_solve(B):
    w = 1
    Y = []
    while len(Y) == 0:
        p = 1
        while p <= w:
            X = enumerate_vectors_weight_p(B.ncols(), p)
            for c_1 in X:
                c_2 = - c_1 * B.T
                if c_1.hamming_weight() + c_2.hamming_weight() <= w:
                    c = vector(Fq, list(c_1) + list(c_2))
                    Y.append(c) #REMEMBER TO REVERSE THE PERMUTATION!
            p += 1
        w += 1
    return Y

def stern_solve(A, l):
    Y = []
    A1 = A[:, :ceil(A.ncols()/2)]
    A2 = A[:, ceil(A.ncols()/2): ]
    w = 1

    while len(Y) == 0:
        p=1
        while p <= w:
            #Building first list
            X1 = enumerate_vectors_weight_p(ceil(A.ncols()/2), p)
            L1 = [ [x1 * A1.T, x1] for x1 in X1] 
            
            #Building second list
            if A.ncols() % 2 == 0:
                X2 = X1
            else:
                X2 = enumerate_vectors_weight_p(floor(A.ncols()/2), w-p)
            L2 = [ [-x2 * A2.T, x2] for x2 in X2] 

            #Sorting both lists
            L1.sort()
            L2.sort()

            #Performing Binary search
            values1 = {tuple(entry[0]) for entry in L1}
            values2 = {tuple(entry[0]) for entry in L2}
            matching_values = values1.intersection(values2)
            for value in matching_values:
                Y.append(vector(Fq, list(L1[binary_search(L1,vector(Fq, value))][1])+list(L2[binary_search(L2, vector(Fq,value))][1])))

            p += 1
        w += 1
    return Y   

def enumerate_vectors_weight_p(size, p):
    v = [0] * (size-p) + [1] * p
    perms = Permutations(v, len(v)).list()
    vectors = [vector(Fq, perm) for perm in perms]
    return vectors

def isd_lee_brickell(H):
    # Cycle until a valid permutation matrix P is found
    while True:
        P = sample_permutation_matrix(n)
        #PH = H*P
        PH = P*H
        if check_random_permutation(PH, 0):
            Pinv = P.T
            break

    PH = pge(PH, 0)
    B = PH[:, :k]
    Y = lee_brickell_solve(B)
    Z = [ y * Pinv for y in Y]

    return Z

def isd_stern(H, l):
    # Cycle until a valid permutation matrix P is found
    while True:
        P = sample_permutation_matrix(n)
        PH = H*P
        if check_random_permutation(PH, l):
            Pinv = P.T
            break

    PH = pge(PH, l)
    A = PH[:l , :k+l]
    B = PH[l:, :k+l]
    Y = stern_solve(A, l)
    Z = [ vector(Fq, list(y) + list(-y*B.T)) * Pinv for y in Y]

    return Z

def binary_search(sorted_list, target_code):
    left, right = 0, len(sorted_list) - 1

    while left <= right:
        mid = (left + right) // 2
        current_code = sorted_list[mid][0]

        if current_code == target_code:
            return mid  # Code found at index mid
        elif current_code < target_code:
            left = mid + 1
        else:
            right = mid - 1

    return None  # Code not found


'''
ALL STEPS (for solving a random code with Lee and Brickell):
    1) [x] Generate the random code (maybe generate randomly a parity check matrix)
    2) [x] Generate a valid random permutation (repeat until it's found)
    3) [x] Perform PGE on the permuted parity check matrix and obtain B (an (n-k)x(k) matrix => (r)x(k) matrix)
    4) [x] Define a low weight p to start itering on
    5) [x] Enumerate all possible c' (which size is a (1)x(k) vector)
    6) [x] For each generated c' compute c''=-c'*B^(T)
    7) [x] If the sum of the weights of c' and c'' is w, add it to the set of solutions
    8) [x] If the set of solutions is not null return it, otherwise increase w of 1 and go back to step 4
    9) [ ] Calculate time complexity

Some random thoughts:
    - what is usually a good starting w (i guess in relation to the code size)
    - what are the typical code sizes
    - what is a good value for p? is it iterated too?

ALL STEPS (for solving a random code with Stern):
    1) [x] Generate the random code (maybe generate randomly a parity check matrix)
    2) [x] Generate a valid random permutation (repeat until it's found)
    3) [x] Perform PGE on the permuted parity check matrix and obtain A and B, thus A1 and A2
    4) [x] Define a low weight p to start itering on
    5) [x] Enumerate all possible x' and x''
    6) [X] Generate lists L1 and L2
    7) [ ] Find matches in a proper way (complexity speaking)
    8) [x] For each match c' calculate c'' =  -c'*B.T
    8) [x] If the set of solutions is not null return it, otherwise increase w of 1 and go back to step 4
    ?) [ ] Use an optimal l
    9) [ ] Calculate time complexity
'''

if __name__ == "__main__":
    print("----- ISD Algorithm -----")
    print("Solve version: Brickell and Lee")
    print("Code properties:")
    print(f"- Length: {n}")
    print(f"- Dimension: {k}")
    print(f"- Redundancy: {r}")
    print("Starting solve...")
    #C = random_linear_code(Fq, n, k)
    #H = C.parity_check_matrix()
    H = generate_random_ldpc(11,22,6,3)
    #minimum_distance_codewords = isd_lee_brickell(H)
    #sample_codeword = minimum_distance_codewords[0]
    #minimum_distance = sample_codeword.hamming_weight()
    #print("DONE!")
    #print(f"Checking solutions...")
    #for codeword in minimum_distance_codewords:
    #    if codeword * H.T != 0:
    #            print("A codeword doesn't doesn't belong to the code!")
    #print(f"Minimum distance: {minimum_distance}")
    #print(f"Minimum distance codewords found: {len(minimum_distance_codewords)}")
    #print(f"First entry: {sample_codeword}")
    #print(f"Is result correct? {C.minimum_distance() == minimum_distance}")
    #print(f"Actual minimum_distance: {C.minimum_distance()}")
    #print(f"Checking solutions...")
    #print("---------------------------------------------")
    #print("Solve version: Stern")
    #print("Code properties:")
    #print(f"- Length: {n}")
    #print(f"- Dimension: {k}")
    #print(f"- Redundancy: {r}")
    #print(f"- l variable: {l}")
    #print("Starting solve...")
    minimum_distance_codewords = isd_stern(H, l)
    sample_codeword = minimum_distance_codewords[0]
    minimum_distance = sample_codeword.hamming_weight()
    print("DONE!")
    print(f"Checking solutions...")
    print(minimum_distance_codewords)
    for codeword in minimum_distance_codewords:
        if codeword * H.T != 0:
                print("A codeword doesn't doesn't belong to the code!")
        else:
            distance = codeword.hamming_weight()
            minimum_distance = distance if distance < minimum_distance else minimum_distance
    print(f"Minimum distance: {minimum_distance}")
    print(f"Minimum distance codewords found: {len(minimum_distance_codewords)}")
    print(f"First entry: {sample_codeword}")
    print(f"Is result correct? {C.minimum_distance() == minimum_distance}")
    print(f"Actual minimum_distance: {C.minimum_distance()}")
