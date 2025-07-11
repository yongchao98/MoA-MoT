import math
from fractions import Fraction

def solve_heat_transfer_problem():
    """
    Solves the heat transfer problem based on the provided image and parameters.
    """
    # Step 1: Define n_0 and m_0 based on the analysis of the plots.
    # n_0 is the plot number for the carbon steel fin. Our analysis shows it's plot 9.
    n_0 = 9
    # m_0 is 1 for a circular fin and -1 for a square fin. Our analysis
    # shows the carbon steel fin is circular.
    m_0 = 1

    print(f"Based on analysis of the plots:")
    print(f"The carbon steel fin is in plot n_0 = {n_0}")
    print(f"The geometry of the carbon steel fin is circular, so m_0 = {m_0}\n")

    # Step 2: Calculate R(c) for the circular geometry.
    # Given conditions: hL/k = ln(13), 4L/d = ln(13)
    # This leads to mL = ln(13) and h/(mk) = 1.
    # R(c) = 1 / tanh(ln(13))
    
    # Using fractions for precision
    # tanh(ln(x)) = (x - 1/x) / (x + 1/x) = (x^2 - 1) / (x^2 + 1)
    x_c = 13
    tanh_ln13_num = x_c**2 - 1
    tanh_ln13_den = x_c**2 + 1
    tanh_ln13 = Fraction(tanh_ln13_num, tanh_ln13_den) # 168/170 = 84/85
    
    R_c = 1 / tanh_ln13 # 85/84
    
    print("For the circular case:")
    print(f"R(c) = 1 / tanh(ln(13)) = 1 / ({tanh_ln13.numerator}/{tanh_ln13.denominator}) = {R_c.numerator}/{R_c.denominator}\n")

    # Step 3: Calculate R(s) for the square geometry.
    # Given conditions: hL/k = ln(2), 4L/w = ln(2)
    # This leads to mL = ln(2) and h/(mk) = 1.
    # R(s) = 1 / tanh(ln(2))

    x_s = 2
    tanh_ln2_num = x_s**2 - 1
    tanh_ln2_den = x_s**2 + 1
    tanh_ln2 = Fraction(tanh_ln2_num, tanh_ln2_den) # 3/5

    R_s = 1 / tanh_ln2 # 5/3

    print("For the square case:")
    print(f"R(s) = 1 / tanh(ln(2)) = 1 / ({tanh_ln2.numerator}/{tanh_ln2.denominator}) = {R_s.numerator}/{R_s.denominator}\n")

    # Step 4: Compute the final expression.
    # Expression: n_0 * (R(c) / R(s))^m_0
    
    # Calculate R(c) / R(s)
    ratio_R = R_c / R_s

    # Final result
    result = n_0 * (ratio_R) ** m_0

    # Print the final equation with all numbers substituted
    print("Final Calculation:")
    final_equation = f"{n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^{m_0}"
    print(f"Expression: n_0 * (R(c) / R(s))^m_0")
    print(f"Substituting values: {final_equation}")
    print(f"= {n_0} * ({ratio_R.numerator}/{ratio_R.denominator})^{m_0}")
    print(f"= {result.numerator}/{result.denominator}")
    
    # Output the final answer in the required format
    print(f"\nFinal Answer: {result.numerator}/{result.denominator}")


solve_heat_transfer_problem()