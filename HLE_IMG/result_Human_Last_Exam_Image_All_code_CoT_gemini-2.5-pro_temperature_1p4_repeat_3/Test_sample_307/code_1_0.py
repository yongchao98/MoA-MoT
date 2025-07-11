import sympy.physics.wigner as wigner
import numpy as np

def get_infinity_norm(j1, j2, j3):
    """
    Calculates the infinity-norm for a Wigner 3-j symbol by finding the
    maximum absolute value over all allowed m-values.
    """
    max_abs_val = 0.0
    
    # Iterate over all possible m values (which are integers for integer j)
    m1_range = np.arange(-j1, j1 + 1)
    m2_range = np.arange(-j2, j2 + 1)
    
    for m1 in m1_range:
        for m2 in m2_range:
            # From the selection rule m1 + m2 + m3 = 0
            m3 = -m1 - m2
            
            # Check the triangle rule for m-values
            if abs(m3) > j3:
                continue
            
            # Calculate the Wigner 3-j symbol
            val = wigner.wigner_3j(j1, j2, j3, m1, m2, m3)
            
            # Update the maximum absolute value found so far
            val_f = abs(float(val))
            if val_f > max_abs_val:
                max_abs_val = val_f
                
    return max_abs_val

def solve_wigner_ratio():
    """
    Solves the problem by calculating norms for a plausible set of j values
    and finding the ratio of the max to the min norm.
    """
    j1 = 8
    j2 = 8
    # The 9 j3 values that result in non-zero symbols for j1=j2=8
    j3_values = [0, 2, 4, 6, 8, 10, 12, 14, 16]

    infinity_norms = []
    for j3 in j3_values:
        norm = get_infinity_norm(j1, j2, j3)
        infinity_norms.append(norm)

    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)
    
    ratio = max_norm / min_norm
    
    print(f"Maximum Infinity-Norm: {max_norm}")
    print(f"Minimum Infinity-Norm: {min_norm}")
    print(f"Ratio of Maximum to Minimum:")
    print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == "__main__":
    solve_wigner_ratio()
