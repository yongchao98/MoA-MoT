import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import N

def calculate_infinity_norm_for_wigner_3j(j1, j2, j3):
    """
    Constructs the Wigner 3-j symbol matrix for a given (j1, j2, j3)
    and calculates its infinity-norm.

    The matrix is indexed by m1 (rows) and m2 (columns).
    """
    # m values range from -j to j
    m1_vals = np.arange(-j1, j1 + 1)
    m2_vals = np.arange(-j2, j2 + 1)
    
    # Initialize the matrix
    matrix = np.zeros((len(m1_vals), len(m2_vals)))
    
    # Populate the matrix with Wigner 3-j symbol values
    for i, m1 in enumerate(m1_vals):
        row_sum = 0.0
        for j, m2 in enumerate(m2_vals):
            # m3 is fixed by the selection rule m1 + m2 + m3 = 0
            m3 = -m1 - m2
            
            # The wigner_3j function from sympy returns 0 if selection rules are not met.
            # We convert the symbolic result from sympy to a floating-point number.
            val = float(N(wigner_3j(j1, j2, j3, m1, m2, m3)))
            matrix[i, j] = val

    # The infinity-norm is the maximum absolute row sum.
    # np.sum(np.abs(matrix), axis=1) calculates the sum of absolute values for each row.
    # np.max() finds the maximum value among these row sums.
    infinity_norm = np.max(np.sum(np.abs(matrix), axis=1))
    
    return infinity_norm

# Based on deduction, the 9 plots correspond to j1=4, j2=4,
# and j3 ranging from 0 to 8.
j1 = 4
j2 = 4
j3_values = range(j1 + j2 + 1)

# Store the calculated norms
norms = []
for j3 in j3_values:
    norm = calculate_infinity_norm_for_wigner_3j(j1, j2, j3)
    norms.append(norm)

# Find the maximum and minimum norms
max_norm = max(norms)
min_norm = min(norms)

# Calculate the ratio
ratio = max_norm / min_norm

# Print the final equation as requested
print("The maximum infinity-norm found is:", max_norm)
print("The minimum infinity-norm found is:", min_norm)
print("\nThe equation for the ratio is:")
print(f"{max_norm} / {min_norm} = {ratio}")
