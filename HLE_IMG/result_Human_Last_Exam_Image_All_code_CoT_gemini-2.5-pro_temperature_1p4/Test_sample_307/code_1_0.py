import numpy as np
from sympy.physics.wigner import wigner_3j

def solve_wigner_ratio():
    """
    This script calculates the ratio of the maximum to the minimum infinity-norm
    among the nine Wigner 3-j symbols shown in the image.
    """
    
    # Based on the image source and visual analysis, the parameters for the 9 plots are:
    # j1 is fixed at 5.
    # j2 corresponds to the row of the plot grid (4, 5, 6).
    # j3 corresponds to the column of the plot grid (4, 5, 6).
    param_list = [
        (5, 4, 4), (5, 4, 5), (5, 4, 6),  # Top row (plots 1, 2, 3)
        (5, 5, 4), (5, 5, 5), (5, 5, 6),  # Middle row (plots 4, 5, 6)
        (5, 6, 4), (5, 6, 5), (5, 6, 6)   # Bottom row (plots 7, 8, 9)
    ]

    infinity_norms = []

    # Iterate through each parameter set to calculate the corresponding infinity-norm
    for params in param_list:
        j1, j2, j3 = params
        
        # Define the ranges for m1 and m2
        m1_min, m1_max = -int(j1), int(j1)
        m2_min, m2_max = -int(j2), int(j2)
        
        num_rows = 2 * m1_max + 1
        num_cols = 2 * m2_max + 1
        
        # Initialize the matrix
        wigner_matrix = np.zeros((num_rows, num_cols))

        # Populate the matrix with Wigner 3-j symbol values
        for i, m1 in enumerate(range(m1_min, m1_max + 1)):
            for k, m2 in enumerate(range(m2_min, m2_max + 1)):
                m3 = -m1 - m2
                # The wigner_3j function returns 0 if selection rules are not met.
                # We can evaluate the symbol directly.
                value = wigner_3j(j1, j2, j3, m1, m2, m3)
                # Convert the sympy object to a float
                wigner_matrix[i, k] = float(value.evalf())
        
        # Calculate the infinity-norm (maximum absolute row sum)
        # ord=np.inf is not implemented for non-square matrices in some older numpy versions.
        # We calculate it manually for robustness.
        inf_norm = np.max(np.sum(np.abs(wigner_matrix), axis=1))
        infinity_norms.append(inf_norm)

    # Find the maximum and minimum norms from the list
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)
    
    # Calculate the ratio
    ratio = max_norm / min_norm
    
    # Print the final result in the requested format
    print(f"Maximum infinity-norm: {max_norm}")
    print(f"Minimum infinity-norm: {min_norm}")
    print(f"Ratio = {max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    solve_wigner_ratio()
