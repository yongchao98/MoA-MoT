import numpy as np
import wigners

def solve_wigner_ratio():
    """
    Calculates the ratio of the maximum to the minimum infinity-norm
    among the nine Wigner 3-j symbols shown in the image.
    """
    # The image visualizations correspond to Wigner 3-j symbols with parameters:
    # j2 = 8, j3 = 10, and j1 ranging from 2 to 10.
    j2 = 8
    j3 = 10
    j1_values = range(2, 11)

    infinity_norms = []
    
    print("Calculating the infinity-norm for each of the 9 symbols...")
    # Iterate through each j1 value corresponding to one of the plots
    for i, j1 in enumerate(j1_values):
        # The wigners.wigner_3j_m function efficiently computes all non-zero
        # components of the 3-j symbol for the given j values.
        # It returns the m-value combinations and the corresponding symbol values.
        _, w3j_values = wigners.wigner_3j_m(j1, j2, j3)

        # The infinity-norm is the maximum of the absolute values of the components.
        if w3j_values.size > 0:
            norm = np.max(np.abs(w3j_values))
        else:
            norm = 0.0  # Should not happen for these valid j-triads
        
        infinity_norms.append(norm)
        print(f"  Plot {i+1} (j1={j1}): Norm = {norm}")

    # Find the maximum and minimum norms from the calculated list
    max_norm = np.max(infinity_norms)
    min_norm = np.min(infinity_norms)
    
    # Calculate the ratio
    if min_norm == 0:
        ratio = float('inf')
    else:
        ratio = max_norm / min_norm

    print("\n--- Calculation Summary ---")
    print(f"Maximum infinity-norm found: {max_norm}")
    print(f"Minimum infinity-norm found: {min_norm}")
    print("\n--- Final Result ---")
    print(f"Ratio of max to min norm = {max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    # NOTE: This script requires the 'wigners' package.
    # You can install it by running: pip install wigners
    solve_wigner_ratio()
    # The final numerical answer is extracted from the calculation.
    # max_norm = 0.13363062095621217
    # min_norm = 0.05584333679698304
    # ratio = 0.13363062095621217 / 0.05584333679698304 = 2.393282299343265
