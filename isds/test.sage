import random
from sage.coding.code_constructions import random_linear_code

Fq = GF(2)
C = random_linear_code(domain, n, k)
return C.parity_check_matrix()
H = Matrix(m,n)
row = random.choice(available_indices)
H[row, i]=1
available_indices.remove(row)
total_ones_indices[row].append(i)
sigma = Permutations(n).random_element()
submatrix = PH[:, -(r-l):]
submatrix.rank()
A[affected_row, i] = A[affected_row][i] + A[addend_row][i]
S = identity_matrix(Fq, r)
P = P[-l:,:].stack(P[:r-l,:])
c_1.hamming_weight() + c_2.hamming_weight() <= w:
c = vector(Fq, list(c_1) + list(c_2))
A1 = A[:, :ceil(A.ncols()/2)]
L2.sort()

values1 = {tuple(entry[0]) for entry in L1}
matching_values = values1.intersection(values2)
Pinv = P.T
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
