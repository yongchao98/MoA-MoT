import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import N

def calculate_wigner_matrix(j1, j2, j3):
    """
    Calculates the matrix of Wigner 3-j symbols for given j values.
    The matrix elements are W(j1, j2, j3, m1, m2, -m1-m2).
    """
    m1_range = np.arange(-j1, j1 + 1)
    m2_range = np.arange(-j2, j2 + 1)
    
    # Initialize a matrix to store the symbol values
    matrix = np.zeros((len(m1_range), len(m2_range)), dtype=float)
    
    # Loop over all possible m1 and m2 values
    for i, m1 in enumerate(m1_range):
        for j, m2 in enumerate(m2_range):
            m3 = -m1 - m2
            # The symbol is non-zero only if |m3| <= j3
            if abs(m3) <= j3:
                # Use sympy to calculate the Wigner 3-j symbol and convert to float
                val = N(wigner_3j(j1, j2, j3, m1, m2, m3))
                matrix[i, j] = float(val)
                
    return matrix

def calculate_infinity_norm(matrix):
    """
    Calculates the infinity-norm of a matrix, which is the maximum
    absolute row sum.
    """
    # Sum the absolute values across columns for each row
    row_sums = np.sum(np.abs(matrix), axis=1)
    # Return the maximum of these sums
    return np.max(row_sums)

def solve_wigner_ratio():
    """
    Main function to solve the problem.
    """
    # Parameters identified from the image source
    j1 = 8
    j2 = 8
    j3_values = range(4, 13) # Corresponds to plots 1 through 9

    norms = []
    
    print("Calculating the infinity-norm for each of the 9 matrices...")
    for j3 in j3_values:
        # 1. Generate the matrix
        wigner_matrix = calculate_wigner_matrix(j1, j2, j3)
        
        # 2. Calculate its infinity-norm
        norm = calculate_infinity_norm(wigner_matrix)
        norms.append(norm)
        # Optional: print the norm for each matrix
        # print(f"  j3={j3}: Norm = {norm:.6f}")

    # 3. Find the maximum and minimum norms from the list
    max_norm = max(norms)
    min_norm = min(norms)

    # 4. Calculate the final ratio
    ratio = max_norm / min_norm

    print("\n--- Final Calculation ---")
    print(f"The maximum infinity-norm found is: {max_norm}")
    print(f"The minimum infinity-norm found is: {min_norm}")
    print(f"The ratio is {max_norm} / {min_norm} = {ratio}")

if __name__ == "__main__":
    solve_wigner_ratio()