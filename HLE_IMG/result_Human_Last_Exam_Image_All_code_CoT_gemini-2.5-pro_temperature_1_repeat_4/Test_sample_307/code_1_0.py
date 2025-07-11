import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import N

def calculate_wigner_3j_infinity_norm(j1, j2, j3):
    """
    Constructs the Wigner 3-j symbol matrix for a given (j1, j2, j3)
    and calculates its infinity-norm.
    The matrix is conceptually indexed by m1 (rows) and m2 (columns).
    """
    # Define the ranges for m1, which corresponds to the rows of the matrix.
    m1_vals = np.arange(-j1, j1 + 1)
    m2_vals = np.arange(-j2, j2 + 1)

    max_abs_row_sum = 0.0

    # Iterate over each row (each value of m1) to find the max absolute row sum.
    for m1 in m1_vals:
        current_abs_row_sum = 0.0
        # Sum the absolute values of the elements in this row (fixed m1).
        for m2 in m2_vals:
            # The third m-quantum number is constrained by m1 + m2 + m3 = 0.
            m3 = -m1 - m2
            
            # The symbol is non-zero only if the selection rule |m3| <= j3 is met.
            if abs(m3) <= j3:
                # Calculate the Wigner 3-j symbol using sympy.
                symbol_val = wigner_3j(j1, j2, j3, m1, m2, m3)
                # Add its absolute value to the current row sum.
                current_abs_row_sum += abs(float(N(symbol_val)))
        
        # Update the overall maximum absolute row sum if the current one is larger.
        if current_abs_row_sum > max_abs_row_sum:
            max_abs_row_sum = current_abs_row_sum
            
    return max_abs_row_sum

def solve_task():
    """
    Main function to solve the problem as described.
    """
    # Parameters for the nine Wigner 3-j symbols visualized in the image.
    j2_val = 10
    j3_val = 10
    j1_vals = range(1, 10)

    # Calculate the infinity-norm for each of the nine matrices.
    infinity_norms = []
    for j1_val in j1_vals:
        norm = calculate_wigner_3j_infinity_norm(j1_val, j2_val, j3_val)
        infinity_norms.append(norm)

    # Find the maximum and minimum infinity-norms.
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)

    # Calculate the final ratio.
    if min_norm == 0:
        ratio = float('inf')
    else:
        ratio = max_norm / min_norm

    # Print the final result, showing each number in the equation.
    print("--- Calculation Results ---")
    print(f"Maximum infinity-norm: {max_norm}")
    print(f"Minimum infinity-norm: {min_norm}")
    print("\n--- Final Equation ---")
    print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    solve_task()