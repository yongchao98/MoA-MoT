import numpy as np

def calculate_writhe(n, o_permutation, x_permutation):
    """
    Calculates the writhe of a grid diagram.
    The grid has 'o's at (i, o_permutation[i-1]) and 'x's at (i, x_permutation[i-1]).
    The permutations are 1-indexed.
    """
    # Inverse of a permutation
    def inverse_permutation(p):
        inv_p = np.zeros(len(p), dtype=int)
        for i, val in enumerate(p):
            inv_p[val-1] = i + 1
        return inv_p

    pi = o_permutation
    sigma = x_permutation
    
    pi_inv = inverse_permutation(pi)
    sigma_inv = inverse_permutation(sigma)
    
    # Calculate vertical and horizontal orientation signs
    v_sign = np.sign(sigma - pi)
    h_sign = np.sign(sigma_inv - pi_inv)

    # Calculate writhe using the summation formula
    # w = (sum of v_signs) * (sum of h_signs) - sum of corrections for marked cells
    total_grid_sum = np.sum(v_sign) * np.sum(h_sign)

    # Correction for 'o' positions (i, pi[i-1])
    sum_o_corr = 0
    for i in range(1, n + 1):
        col_idx = i - 1
        row_idx = pi[col_idx] - 1
        sum_o_corr += v_sign[col_idx] * h_sign[row_idx]

    # Correction for 'x' positions (i, sigma[i-1])
    sum_x_corr = 0
    for i in range(1, n + 1):
        col_idx = i - 1
        row_idx = sigma[col_idx] - 1
        sum_x_corr += v_sign[col_idx] * h_sign[row_idx]

    writhe = total_grid_sum - sum_o_corr - sum_x_corr
    return int(writhe)

# The knot specified by the 5x5 grid is a trefoil knot.
# To find its maximal Thurston-Bennequin number, we must analyze its
# minimal grid diagram, which is a 3x3 grid.

# Grid number for the minimal trefoil diagram
n_min = 3

# For a standard minimal grid diagram of a (right-handed) trefoil, the 'o's are
# on the main diagonal and the 'x's follow a cyclic permutation.
# O positions: (1,1), (2,2), (3,3) -> identity permutation
# X positions: (1,2), (2,3), (3,1) -> cyclic permutation (2,3,1)

o_perm_min = np.array([1, 2, 3])
x_perm_min = np.array([2, 3, 1])

# Calculate the writhe of this minimal 3x3 diagram
w_min = calculate_writhe(n_min, o_perm_min, x_perm_min)

# The maximal Thurston-Bennequin number is tb = w - n for the minimal grid
tb_max = w_min - n_min

# Print the final equation with the computed values
print(f"The maximal Thurston-Bennequin number is calculated from the minimal grid diagram (n={n_min}).")
print("The writhe (w) of this minimal diagram has been calculated.")
print(f"The final calculation is: w - n = {w_min} - {n_min} = {tb_max}")
print("\nFinal Answer:")
print(tb_max)
