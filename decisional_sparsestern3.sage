import random
import logging
import numpy as np
from itertools import combinations

#  Generate LDPC Code
#  Apply right and left permutations

logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)

def generate_matrix():
	M = Matrix(GF(2), r, n, lambda i, j: random.choices([0, 1], weights=[0.95,0.05], k=1)[0])
	while any(col.is_zero() for col in M.columns()):
		M = Matrix(GF(2), r, n, lambda i, j: random.choices([0, 1], weights=[0.80,0.2], k=1)[0])
	return M

def move_l_rows_up_supports_left(H,l):
	support = [0 for _ in range(H.ncols())]
	support_length = 0
	l_rand_rows = random.sample(range(H.nrows()), l)
	logging.debug(f"Random rows: {l_rand_rows}")
	for row in l_rand_rows:
		for col in range(H.ncols()):
			if H[row,col] != 0:
				if support[col] == 0:
				    support[col] = 1
				    support_length += 1

	P = identity_matrix(GF(2),H.ncols())
	col_swap_array = support.copy()
	logging.debug(f"Support length: {support_length}")
	for col in range(support_length,H.ncols()):
		if col_swap_array[col] == 1:
			for swap_col in range(len(col_swap_array)):
				if col_swap_array[swap_col] == 0:
				    # Modify P to swap (col,swap_col)
				    P[col,col] = 0
				    P[col,swap_col]=1
				    P[swap_col,swap_col] = 0
				    P[swap_col,col] = 1
				    col_swap_array[swap_col] = 2; break

	Q = identity_matrix(GF(2),H.nrows())
	row_swap_array = [0 for _ in range(l)]
	for row in l_rand_rows:
		if row < l:
			row_swap_array[row] = 1
	logging.debug(f"Swappable rows: {row_swap_array}")
	for row in l_rand_rows:
		if row >= l:
			for swap_row in range(l):
				if row_swap_array[swap_row] == 0:
					Q[row, row] = 0
					Q[row, swap_row] = 1
					Q[swap_row, swap_row] = 0
					Q[swap_row, row] = 1
					row_swap_array[swap_row] = 2; break
	logging.debug(f"Swappable rows after swap: {row_swap_array}")
	return Q*H*P,P,Q, support_length;

# EVITABILE FACENDOLO DIRETTAMENTE NEL PASSAGGIO PRIMA!
def shift_column_right(H,offset):
	P = identity_matrix(GF(2),H.ncols())
	logging.debug(f"Shifting cols {offset} right")
	for i in range(z):
		P[i,i] = 0
		P[i,i+offset] = 1
	for i in range(offset):
		P[z+i,z+i] = 0
		P[z+i,i] = 1
	return H*P, P

# Note: LDPC have no null columns, we consider pge to be possible only with row sums and permutation
# PGE to create a x size I matrix on bottom right corner
def pge(A,x):
	if A.rank() < x:
		print("Rank 2 low!")
		print(A)
		print(z)
		print(x)
		return
	else:
		print("Rank OK!")

	logging.debug(f"Searching for {x} indip. columns")
	selected_cols = [0]
	current_rank = 1
	for i in range(1,A.ncols()):
		submatrix = A[:,(selected_cols + [i])]
		if submatrix.rank() > current_rank:
			selected_cols.append(i)
			current_rank +=1
			if current_rank == x:
				break
	logging.debug(f"Columns found: {selected_cols}")

	logging.debug(f"Putting columns to rightmost position")
	RP = identity_matrix(GF(2), A.ncols())
	col_swap_array = [0 for _ in range(x)]
	offset = A.ncols()-x
	for col in selected_cols:
		if col >= offset:
			col_swap_array[col-offset] = 1
	logging.debug(f"Swappable cols: {col_swap_array}")
	for col in selected_cols:
		if col < offset:
			for swap_col in range(offset, A.ncols()):
				if col_swap_array[swap_col-offset] == 0:
					RP[col, col] = 0
					RP[col, swap_col] = 1
					RP[swap_col, swap_col] = 0
					RP[swap_col, col] = 1
					col_swap_array[swap_col-offset] = 2; break
	A1=A*RP
	if A1[:,-x:].rank() == x:
		print("Success!")
	else:
		print("Failure!")
		print(A)
		print("====")
		print(RP)
		return

	LP = identity_matrix(GF(2),A1.nrows())
	for i in range(x):
		LP1 = identity_matrix(GF(2),A1.nrows())
		LP2 = identity_matrix(GF(2),A1.nrows())
		pivot_row = A1.nrows()-x+i
		pivot_col = A1.ncols()-x+i
		pivot = A1[pivot_row,pivot_col]
		# Make sure all pivots are 1
		if pivot == 0:
			addable_rows = [j for j in range(0,A1.nrows()-x)]
			addable_rows += [j for j in range(A1.nrows()-x+i,A1.nrows())]
			for row in addable_rows:
				if A1[row,pivot_col] == 1:
					LP1[pivot_row,row] = 1
					A1=LP1*A1
					break
		
		pivot = A1[pivot_row,pivot_col]
		if pivot != 1:
			logging.debug(f"{pivot_row}-{pivot_col}: {pivot}")
			print(A1)
			exit(-1)
		for row in range(A1.nrows()):
			if row == pivot_row:
				continue
			if A1[row,pivot_col] == 1:
				LP2[row,pivot_row] = 1
		A1=LP2*A1
		LP=LP2*LP1*LP

	if A1[-(r-l):,-(r-l):] != identity_matrix(GF(2), r-l):
		print(f"Error creating {r-l} Identity matrix")
		print("=========PGE===========")
		print(A1)
		print("=========PGE===========")
		exit(-1)
	
	return LP,RP

def enum_vec_w_p(size, p):
	X = []
	for comb in combinations(range(size),p):
		x = vector(GF(2),size)
		for j in comb:
			x[j] = 1
		X.append(x)
	return X

def sparseStern(H2):
	A = H2[:l1, :z]
	z1 = floor(z/2)
	z2 = ceil(z/2)
	A1 = A[:,:z1]
	A2 = A[:,z1:]
	B = H2[l1:l, :k+l]
	B1 = B[:, :z]
	B2 = B[:, z:]
	C = H2[-(r-l):, :k+l]
	
	p1 = ceil(z1/25)
	p2 = ceil((k+l-z)/25)
	p3 = p1 + p2
	p = p1 + p2 + p3
	combs11 = combinations(range(z1),p1)
	combs12 = combinations(range(z2),p1)
	combs22 = combinations(range(k+l-z),p2)

	collisions = 0
	lowest_weight = n
	
	logging.debug("Starting collision search")
	for comb11 in combs11:
		x11 = vector(GF(2),z1)
		for j in comb11:
			x11[j] = 1
		for comb12 in combs12:
			x12 = vector(GF(2),z2)
			for j in comb12:
				x12[j] = 1
			if x11*A1.T == -x12*A2.T:
				logging.debug("Found collision list 1")
				for comb2 in combs22:
					x21 = vector(GF(2), list(x11) + list(x12))
					x22 = vector(GF(2), k+l-z)
					for j in comb2:
						x22[j] = 1
					#print(x22)
					#print(x11)
					#print(x12)
					#print(x21)
					if x21*B1.T == -x22*B2.T:
						logging.debug("Found collision list 2")
						collisions += 1;
						x31 = vector(GF(2), list(x21) + list(x22))
						x32 = -x31 * C.T
						x = vector(GF(2), list(x31) + list(x32))
						w = x.hamming_weight()
						lowest_weight = min(w, lowest_weight)
						logging.debug(f"Collisions: {collisions}, searching {p}, lowest {lowest_weight}")
						if w <= p:
							return x
	

l1 = 10
l2 = 20
l = l1 + l2
rate = 2/3
n = 1000
k = ceil(n * rate)
r = n-k
#p1=10
#p2=5

# GENERATE MATRIX
M = generate_matrix()

LEFT = identity_matrix(GF(2),r)
RIGHT = identity_matrix(GF(2),n)

# CREATE H1
QMP,P,Q,z = move_l_rows_up_supports_left(M,l1)
if z > k+l:
	print("ERROR: z > u")
	exit()

LEFT = Q*LEFT
RIGHT = RIGHT*P



# CREATE H2 (PGE ON LAST r-l1 ROWS)
D = QMP[-(r-l1):, -(n-z):]
LP,RP = pge(D, r-l)
Il1 = identity_matrix(GF(2),l1)
Iz = identity_matrix(GF(2),z)

LS = Il1.augment(Matrix(l1,r-l1))
LS = LS.stack(Matrix(r-l1,l1).augment(LP))

RS = Iz.augment(Matrix(z,n-z))
RS = RS.stack(Matrix(n-z,z).augment(RP))

LEFT = LS*LEFT
RIGHT = RIGHT*RS

H1=LS*QMP*RS
logging.debug(f"u-z={k+l-z}") 
if (k+l-z <= 0):
	print('Error: u-z <= 0')
	exit()

#At this point, MEET IN THE MIDDLE!

c = sparseStern(H1)

if c == None:
	print("No codeword found")
else:
	Rinv = RIGHT.inverse()
	Linv = LEFT.inverse()
	print(f"Codeword: {c}")
	print(f"Weight: {c.hamming_weight()}")
	print(f"Parity: {c*Rinv*M.T}")
	print(f"z: {z}, u: {k+l-z}, l1: {l1}, l2: {l2}")
