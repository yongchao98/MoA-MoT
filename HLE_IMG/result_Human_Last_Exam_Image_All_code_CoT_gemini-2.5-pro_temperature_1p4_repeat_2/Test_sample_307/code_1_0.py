import numpy as np
from sympy.physics.wigner import wigner_3j
import sys

def calculate_infinity_norm(j1, j2, j3):
    """
    Calculates the infinity-norm of a Wigner 3-j symbol matrix.
    The matrix is defined by W_m1m2 = (j1 j2 j3 | m1 m2 -m1-m2).
    The infinity-norm is the maximum absolute row sum.
    """
    m1_vals = np.arange(-j1, j1 + 1)
    m2_vals = np.arange(-j2, j2 + 1)
    
    max_row_sum = 0
    
    # Iterate over rows (indexed by m1)
    for m1 in m1_vals:
        current_row_sum = 0
        # Sum absolute values across columns (indexed by m2)
        for m2 in m2_vals:
            m3 = -m1 - m2
            
            # The 3j-symbol is non-zero only if the triangular inequalities and
            # magnetic quantum number rules are satisfied.
            # sympy's wigner_3j handles this automatically, returning 0 if invalid.
            # Adding an explicit check for |m3| <= j3 can speed up the calculation.
            if abs(m3) <= j3:
                val = wigner_3j(j1, j2, j3, m1, m2, m3)
                # The result from sympy can be symbolic, so we evaluate it to a float
                current_row_sum += abs(val.evalf())
                
        if current_row_sum > max_row_sum:
            max_row_sum = current_row_sum
            
    return max_row_sum

def solve_wigner_ratio():
    """
    Main function to solve the problem.
    """
    # As identified, the parameters for the plots are:
    j1 = 8
    j2 = 10
    # Plots 1 to 9 correspond to j3 = 2, 3, ..., 10
    j3_values = range(2, 11)
    
    infinity_norms = []
    
    print("Calculating the infinity-norm for each of the 9 symbols...")
    # This calculation can take a few moments.
    # To provide progress feedback, we print a dot for each completed calculation.
    for j3 in j3_values:
        norm = calculate_infinity_norm(j1, j2, j3)
        infinity_norms.append(norm)
        sys.stdout.write('.')
        sys.stdout.flush()
    
    print("\nCalculation complete.")

    # Find the maximum and minimum norms from the list
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)

    # Calculate the required ratio
    ratio = max_norm / min_norm
    
    print(f"\nMaximum infinity-norm found: {max_norm}")
    print(f"Minimum infinity-norm found: {min_norm}")
    
    print("\nThe ratio is calculated as follows:")
    print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    solve_wigner_ratio()
