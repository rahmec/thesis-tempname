import random
import logging
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
    if size < p:
    	raise ValueError("Size is less than p")
    v = [0] * (size-p) + [1] * p
    perms = Permutations(v, len(v)).list()
    vectors = [vector(GF(2), perm) for perm in perms]
    return vectors

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

def sparseStern(H2):
	A = H2[:l1, k+l-z:k+l]
	z1 = floor(z/2)
	z2 = ceil(z/2)
	A1 = A[:,:z1]
	A2 = A[:,z1:]
	#print("=====A=====")
	#print(A)
	#print("=====A1=====")
	#print(A1)
	#print("=====A2=====")
	#print(A2)
	
	p1 = ceil(z1/10)
	logging.debug(f"Enumerating {z1} pick {p1}")
	X11 = enum_vec_w_p(z1,p1)
	logging.debug("Enumerated x11")
	if z % 2 == 0:
		X12 = X11.copy()
	
	else:
		X12 = enum_vec_w_p(z2,p1)
	logging.debug("Enumerated x12")
	
	logging.debug("Calculating lists")
	L11 = [[x11*A1.T,x11] for x11 in X11]
	L12 = [[-x12*A2.T,x12] for x12 in X12]


	logging.debug("Sorting lists")
	L11.sort()
	L12.sort()

	logging.debug("Searching collisions")
	#Performing First Binary search
	values1 = {tuple(entry[0]) for entry in L11}
	values2 = {tuple(entry[0]) for entry in L12}
	matching_values = values1.intersection(values2)
	logging.debug("Found matches")
	Y1 = []
	for value in matching_values:
		Y1.append(vector(GF(2), 
		list(L11[binary_search(L11,vector(GF(2), value))][1])
		+list(L12[binary_search(L12, vector(GF(2),value))][1])))
	
	logging.debug("First small instance done")
	
	B = H2[l1:l, :k+l]
	B1 = B[:, :(k+l-z)]
	B2 = B[:, k+l-z:]
	#print("=====B=====")
	#print(B)
	#print("=====B1=====")
	#print(B1)
	#print("=====B2=====")
	#print(B2)
	
	p2 = ceil((k+l-z)/10)
	p = p1 + p2
	logging.debug(f"Enumerating {z2} pick {p2}")
	X21 = enum_vec_w_p((k+l-z),p2)
	logging.debug("Enumerataed X21")
	X22 = Y1

	logging.debug("Calculating lists")
	L21 = [[x21*B1.T,x21] for x21 in X21]
	L22 = [[-x22*B2.T,x22] for x22 in X22]
	
	logging.debug("Sorting lists")
	L21.sort()
	L22.sort()

	#Performing Second Binary search
	logging.debug("Searching collisions")
	values1 = {tuple(entry[0]) for entry in L21}
	values2 = {tuple(entry[0]) for entry in L22}
	matching_values = values1.intersection(values2)
	logging.debug("Found matches")
	Y2 = []
	for value in matching_values:
		Y2.append(vector(GF(2), 
		list(L21[binary_search(L21,vector(GF(2), value))][1])
		+list(L22[binary_search(L22, vector(GF(2),value))][1])))

	logging.debug("Second small instance done")
	
	C = H2[-(r-l):, :k+l]
	Y3 = []
	for value in Y2:
		x3 = -value*C.T
		Y3.append(vector(GF(2),list(value)+list(x3)))

	#logging.debug("Rechecking all of them")
	#for x in Y3:
	#	x11 = x[k+l-z:k+l-z2]
	#	x12 = x[k+l-z2:k+l]
	#	x21 = x[:k+l-z]
	#	x22 = x[k+l-z:k+l]
	#	x31 = x[:k+l]
	#	x32 = x[k+l:]
	#	if x11 * A1.T != -x12 * A2.T:
	#		logging.debug("First equation not respected!")
	#	if x21 * B1.T != -x22 * B2.T:
	#		logging.debug("Second equation not respected!")
	#	if x31 * C.T != x32:
	#		logging.debug("Third equation not respected!")
	#	new1 = matrix(GF(2),l1,k+l-z).augment(A1).augment(A2)
	#	new1 = new1.augment(matrix(GF(2),l1,n-k-l))
	#	new2 = B1.augment(B2).augment(matrix(GF(2),l2,n-k-l))
	#	new3 = C.augment(identity_matrix(GF(2),n-k-l))
	#	new = new1.stack(new2).stack(new3)
	#	if new != H2:
	#		logging.debug("Matrixes don't match")
	#		print("==========H2==========")
	#		print(H2)
	#		print("==========new=========")
	#		print(new)
	#		print("======================")
	#		print(f"r-l={r-l}")
	#		exit(-1)
		
		

	logging.debug("SparseStern done")

	
	return Y3
	

l1 = 2
l2 = 3
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
H2, P2 = shift_column_right(H1,k+l-z)
logging.debug(f"u-z={k+l-z}") 
if (k+l-z <= 0):
	print('Error: u-z <= 0')
	exit()

RIGHT = RIGHT * P2

#At this point, MEET IN THE MIDDLE!

codewords = sparseStern(H2)
print("=======")
print(f"Solution: {len(codewords)} codewords")
Rinv = RIGHT.inverse()
Linv = LEFT.inverse()
print(M == Linv * H2 * Rinv)
for x in codewords:
	#print(x)
	print(x*Rinv*M.T)

print("=============")
#print(H2)
