import numpy as np
from sympy.physics.wigner import wigner_3j
import sympy

def calculate_3j_matrix(j1, j2, j3):
    """
    Calculates the matrix of Wigner 3-j symbols for given j1, j2, j3,
    where the matrix is indexed by m1 and m2.
    
    The script requires the 'sympy' and 'numpy' libraries.
    You can install them using: pip install sympy numpy
    """
    # The range of m quantum numbers for a given j is from -j to +j.
    m1_vals = np.arange(-j1, j1 + 1)
    m2_vals = np.arange(-j2, j2 + 1)
    
    # Initialize the matrix to store the 3-j symbol values.
    # The size is (2*j1 + 1) x (2*j2 + 1).
    matrix = np.zeros((len(m1_vals), len(m2_vals)))
    
    # Iterate over all possible m1 and m2 values.
    for i, m1 in enumerate(m1_vals):
        for j, m2 in enumerate(m2_vals):
            # The m3 quantum number is constrained by m1 + m2 + m3 = 0.
            m3 = -m1 - m2
            
            # The sympy function calculates the precise value of the symbol.
            # It correctly returns 0 if the parameters are invalid.
            val = wigner_3j(j1, j2, j3, m1, m2, m3)
            
            # sympy's wigner_3j can return a symbolic expression.
            # We evaluate it to a floating-point number.
            if isinstance(val, sympy.Expr):
                matrix[i, j] = float(val.evalf())
            else:
                matrix[i, j] = float(val)
                
    return matrix

def calculate_infinity_norm(matrix):
    """
    Calculates the infinity-norm of a matrix, which is the maximum
    absolute row sum. Norm_inf(A) = max_i sum_j |A_ij|.
    """
    # Calculate the sum of absolute values for each row (axis=1).
    row_sums = np.sum(np.abs(matrix), axis=1)
    # The infinity-norm is the maximum of these row sums.
    return np.max(row_sums)

# --- Main execution ---

# Based on visual inspection of the plots, we deduce the parameters.
# The matrices are 11x11, so 2*j+1 = 11 -> j=5. We assume j1=j2=5.
# The nine plots correspond to varying j3 from 1 to 9.
j1 = 5
j2 = 5
j3_values = range(1, 10)

infinity_norms = []

# Loop through each of the nine cases (j3 = 1 to 9).
for j3 in j3_values:
    # 1. Calculate the matrix of Wigner 3-j symbols.
    wigner_matrix = calculate_3j_matrix(j1, j2, j3)
    
    # 2. Calculate its infinity-norm.
    norm = calculate_infinity_norm(wigner_matrix)
    infinity_norms.append(norm)

# 3. Find the maximum and minimum norms from the list of nine norms.
max_norm = max(infinity_norms)
min_norm = min(infinity_norms)

# 4. Calculate the required ratio.
ratio = max_norm / min_norm

# 5. Print the final equation with all the numbers, as requested.
print(f"The maximum infinity-norm found is: {max_norm}")
print(f"The minimum infinity-norm found is: {min_norm}")
print(f"The ratio is {max_norm} / {min_norm} = {ratio}")