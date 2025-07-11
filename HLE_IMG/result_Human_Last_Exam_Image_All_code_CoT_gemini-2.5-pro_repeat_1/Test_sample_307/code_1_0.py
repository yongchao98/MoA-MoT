import sympy
from sympy.physics.wigner import wigner_3j
import numpy as np

def calculate_wigner_infinity_norm(j1, j2, j3):
    """
    Calculates the infinity-norm of a Wigner 3-j symbol matrix.
    The infinity-norm is the maximum absolute value of any element.
    """
    # Define the ranges for the magnetic quantum numbers m1 and m2
    m1_range = range(-j1, j1 + 1)
    m2_range = range(-j2, j2 + 1)
    
    max_abs_val = 0.0
    
    # Iterate through all possible m1 and m2 values
    for m1 in m1_range:
        for m2 in m2_range:
            # The third magnetic quantum number is fixed by the selection rule m1+m2+m3=0
            m3 = -m1 - m2
            
            # The symbol is non-zero only if |m3| <= j3
            if abs(m3) <= j3:
                # Calculate the Wigner 3-j symbol using sympy
                val = wigner_3j(j1, j2, j3, m1, m2, m3)
                
                # Convert the symbolic result to a floating-point number
                float_val = float(val.evalf())
                
                # Update the maximum absolute value found so far
                if abs(float_val) > max_abs_val:
                    max_abs_val = abs(float_val)
                    
    return max_abs_val

def solve_ratio():
    """
    Finds the ratio of the maximum to the minimum infinity-norm
    for the nine Wigner 3-j symbols shown in the image.
    """
    # Parameters identified from the source of the visualization
    j1 = 6
    j2 = 6
    # The nine plots correspond to j3 = 0, 1, 2, ..., 8
    j3_values = range(9)
    
    infinity_norms = []
    print("Calculating infinity-norm for each symbol...")
    for j3 in j3_values:
        norm = calculate_wigner_infinity_norm(j1, j2, j3)
        infinity_norms.append(norm)
        print(f"Symbol with j3={j3} (Plot {j3+1}): Infinity-norm = {norm}")
        
    # Find the maximum and minimum norms from the calculated list
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)
    
    # Calculate the ratio
    ratio = max_norm / min_norm
    
    # Print the final results, including the equation
    print("\n--- Results ---")
    print(f"Maximum infinity-norm: {max_norm}")
    print(f"Minimum infinity-norm: {min_norm}")
    print(f"Ratio = Maximum Norm / Minimum Norm")
    print(f"Ratio = {max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    solve_ratio()
