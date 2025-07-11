import numpy as np

def get_rank(matrix):
    """Computes the rank of a matrix over Z_2."""
    # Perform Gaussian elimination over Z_2
    m, n = matrix.shape
    rank = 0
    pivot_row = 0
    for j in range(n): # Iterate through columns
        if pivot_row < m:
            i = pivot_row
            while i < m and matrix[i, j] == 0:
                i += 1
            
            if i < m:
                # Swap rows i and pivot_row
                matrix[[i, pivot_row]] = matrix[[pivot_row, i]]
                
                # Eliminate other 1s in the same column
                for k in range(m):
                    if k != pivot_row and matrix[k, j] == 1:
                        matrix[k, :] = (matrix[k, :] + matrix[pivot_row, :]) % 2
                pivot_row += 1
    return pivot_row

# For the Klein Bottle, we can use a simple cell complex with
# C_2: 1 face (f)
# C_1: 2 edges (a, b)
# C_0: 1 vertex (v)

# The chain groups over Z_2 are C_2 = Z_2, C_1 = (Z_2)^2, C_0 = Z_2.
dim_c2 = 1
dim_c1 = 2
dim_c0 = 1

# Boundary map d_2: C_2 -> C_1
# The boundary of the face is a+b-a+b = 2b, which is 0 in Z_2.
# So, d_2(f) = 0.
# The matrix for d_2 is a 2x1 zero matrix.
d2_matrix = np.array([[0], [0]])

# Boundary map d_1: C_1 -> C_0
# The boundaries of the edges are d_1(a) = v-v = 0 and d_1(b) = v-v = 0.
# The matrix for d_1 is a 1x2 zero matrix.
d1_matrix = np.array([[0, 0]])

# H_1 = Ker(d_1) / Im(d_2)
# dim(H_1) = dim(Ker(d_1)) - dim(Im(d_2))

# dim(Ker(d_1)) = dim(C_1) - rank(d_1)
rank_d1 = get_rank(d1_matrix.copy())
dim_ker_d1 = dim_c1 - rank_d1

# dim(Im(d_2)) = rank(d_2)
rank_d2 = get_rank(d2_matrix.copy())
dim_im_d2 = rank_d2

# Dimension of the first homology group H_1
dim_h1 = dim_ker_d1 - dim_im_d2

# The group H_1 is (Z_2)^dim_h1, so its size is 2^dim_h1.
size_h1 = 2**dim_h1

# Print the final equation
base = 2
exponent = dim_h1
result = size_h1
print(f"{base}^{exponent} = {result}")
