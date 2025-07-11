import sympy
from sympy.physics.wigner import wigner_3j
import numpy as np

def calculate_3j_norms():
    """
    Calculates the infinity-norms for a family of Wigner 3-j symbols
    and finds the ratio of the maximum to the minimum norm.
    """
    j1 = 5.0
    j2 = 4.0
    
    # The 9 symbols correspond to j3 from |j1-j2| to j1+j2
    j3_values = range(int(abs(j1 - j2)), int(j1 + j2) + 1)
    
    m1_range = np.arange(-j1, j1 + 1)
    m2_range = np.arange(-j2, j2 + 1)
    
    inf_norms = []
    
    print("Calculating infinity-norms for Wigner 3-j symbols (j1={}, j2={}, j3=1..9):\n".format(j1, j2))

    for j3 in j3_values:
        max_abs_w = 0
        for m1 in m1_range:
            for m2 in m2_range:
                m3 = -m1 - m2
                # The wigner_3j function handles selection rules, returning 0 if they are not met.
                # Specifically, it is non-zero only if |m3| <= j3.
                symbol = wigner_3j(j1, j2, j3, m1, m2, m3)
                if symbol != 0:
                    val = abs(float(symbol.evalf()))
                    if val > max_abs_w:
                        max_abs_w = val
        
        inf_norms.append(max_abs_w)
        print("For j3 = {}, the infinity-norm is {:.6f}".format(j3, max_abs_w))
        
    # Find the maximum and minimum norms from the calculated list
    max_norm = max(inf_norms)
    min_norm = min(inf_norms)
    
    # The maximum norm corresponds to j3=1, the minimum to j3=9
    max_norm_j3 = j3_values[inf_norms.index(max_norm)]
    min_norm_j3 = j3_values[inf_norms.index(min_norm)]

    print("\n--------------------------------------------------")
    print("Maximum norm occurs at j3={}.".format(max_norm_j3))
    print("Minimum norm occurs at j3={}.".format(min_norm_j3))
    print("--------------------------------------------------\n")

    # Calculate the ratio
    ratio = max_norm / min_norm
    
    print("The final calculation is:")
    print("Ratio = Maximum Norm / Minimum Norm")
    print("Ratio = {} / {}".format(max_norm, min_norm))
    print("Ratio = {}".format(ratio))

calculate_3j_norms()