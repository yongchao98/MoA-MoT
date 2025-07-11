import math

def solve():
    """
    This function prints the formula for the probability P_m.
    """
    
    # Based on the analysis, the number of favorable pairs (i, j) is N_m = 2*m^2 + 1.
    # The total number of pairs (i, j) is T_m = C(4m+2, 2) = (2m+1)(4m+1) = 8*m^2 + 6*m + 1.
    # The probability P_m = N_m / T_m.
    
    numerator_coeff_m_sq = 2
    numerator_coeff_m = 0
    numerator_const = 1
    
    denominator_coeff_m_sq = 8
    denominator_coeff_m = 6
    denominator_const = 1

    # The problem asks to find the formula for P_m.
    # We will print the formula as a string.
    
    numerator_str = f"({numerator_coeff_m_sq}*m**2 + {numerator_const})"
    denominator_str = f"({denominator_coeff_m_sq}*m**2 + {denominator_coeff_m}*m + {denominator_const})"
    
    print("The formula for the probability P_m is:")
    print(f"P_m = {numerator_str} / {denominator_str}")

solve()