import sympy
from sympy.physics.wigner import wigner_3j
import numpy as np

def solve_wigner_ratio():
    """
    Calculates the ratio of the maximum to minimum infinity-norm for the nine
    Wigner 3-j symbols specified by the problem.
    """
    j1 = 8
    j2 = 8
    m3 = 0
    infinity_norms = []
    
    # The infinity-norm for the matrix A^l is max_{m1} sum_{m2} |(8 8 l; m1 m2 0)|.
    # Due to the m1+m2+m3=0 selection rule, this simplifies to max_{m1} |(8 8 l; m1 -m1 0)|.

    print("Calculating the infinity-norm for each of the 9 plots:")
    for l in range(1, 10):
        j3 = l
        max_abs_for_l = 0
        # Iterate through all possible m1 values (-j1 to j1)
        for m1 in range(-j1, j1 + 1):
            m2 = -m1
            # Calculate the Wigner 3-j symbol using sympy
            # The result is a sympy object, so we evaluate it to a float with .n()
            val = wigner_3j(j1, j2, j3, m1, m2, m3).n()
            abs_val = abs(val)
            if abs_val > max_abs_for_l:
                max_abs_for_l = abs_val
        
        infinity_norms.append(max_abs_for_l)
        print(f"Plot {l} (j3={l}): Norm = {max_abs_for_l:.6f}")

    # Find the maximum and minimum norms from the list
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)
    
    max_norm_l = infinity_norms.index(max_norm) + 1
    min_norm_l = infinity_norms.index(min_norm) + 1
    
    print(f"\nMaximum norm found: {max_norm:.6f} (from plot {max_norm_l})")
    print(f"Minimum norm found: {min_norm:.6f} (from plot {min_norm_l})")

    # Calculate the final ratio
    ratio = max_norm / min_norm
    
    print("\nThe ratio of the maximum to the minimum infinity-norm is:")
    print(f"{max_norm:.6f} / {min_norm:.6f} = {ratio:.4f}")

solve_wigner_ratio()