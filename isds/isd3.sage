import random
from sage.coding.code_constructions import random_linear_code

# Remember to adjust notation to the inernational pyhon standard before goin onward

# Define the various sizes

p0 = 1
n = 22
r = 11
k = n - r
l = 5
Fq = GF(2)

def bin_entropy(x):
    
    return -x*log(x,2)-(1-x)*log(1-x,2)

#equation to find the value of t, as a function of theta_val, ell and alpha_val
def find_t_numerical(theta_val, ell, alpha_val):

    R.<x> = RR[]
    f = (1-2*x)*log((1-x)/x)-2 - 2*theta_val*ell/alpha_val #function whose root returns t
    
    t = find_root(f, 1e-6, 0.5-1e-6) #numerical root finding
    
    return t

#compute min distance
#n, alpha=1-R, ell=colweight
def min_distance_ldpc(n, alpha_val, ell):
    
    b_theta = -1
    w = 2
    while b_theta <0:
        w += 1
        theta_val = w/n
        t = find_t_numerical(theta_val, ell, alpha_val)
        
        p_theta = -alpha_val*log(2)+alpha_val*bin_entropy(t)+theta_val*ell*log(1-2*t)
        
        b_theta = bin_entropy(theta_val) + p_theta
        
    return w

def generate_random_linear_code(domain, n, k):
    C = random_linear_code(domain, n, k)
    return C.parity_check_matrix()

def generate_random_ldpc(r,n,wr,wc):
    return generate_matrix(r,n,wr,wc)

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

if __name__ == "__main__":
    print("----- ISD Algorithm -----")
    # Iterate over r, n, wr, wc, l
    for n in rang
    H = generate_random_ldpc(r,n,6,3)
    minimum_distance_codewords = isd_stern(H, l)
    sample_codeword = minimum_distance_codewords[0]
    minimum_distance = sample_codeword.hamming_weight()
