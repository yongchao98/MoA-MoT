import math
from fractions import Fraction

def solve_fin_problem():
    """
    Solves the heat transfer fin problem based on the provided image and text.
    """

    # Step 1 & 2: Calculate R(c) for the circular fin
    # tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85
    tanh_ln13 = Fraction(84, 85)
    # R(c) = 1 / tanh(ln(13))
    R_c = 1 / tanh_ln13
    
    # Step 3: Calculate R(s) for the square fin
    # tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5
    tanh_ln2 = Fraction(3, 5)
    # R(s) = 1 / tanh(ln(2))
    R_s = 1 / tanh_ln2

    # Step 4: Identify n_0 and m_0 from the plots
    # Based on the detailed analysis in the plan:
    # - The Q_circ/Q_sq ratio is consistent for Copper (~1.095) and one pairing of Platinum fins (Q(5)/Q(3) ~ 1.112).
    # - This consistency implies that Plot 9 represents the Carbon Steel fin.
    # - Therefore, n_0 = 9.
    # - The geometry of the fin in Plot 9 is identified as square by pairing it with the 3D model in Plot 4.
    # - For a square geometry, m_0 = -1.
    n_0 = 9
    m_0 = -1
    
    # Step 5: Calculate the final expression: n_0 * (R(c)/R(s))^m_0
    ratio_R = R_c / R_s
    final_value = n_0 * (ratio_R)**m_0

    # Output the components of the final equation and the result
    print("This script calculates the value of n_0 * (R(c)/R(s))^m_0")
    print("\n--- Intermediate Values ---")
    print(f"n_0 (plot number of carbon steel fin) = {n_0}")
    print(f"m_0 (geometry parameter for carbon steel fin) = {m_0}")
    print(f"R(c) = 1 / tanh(ln(13)) = {R_c.numerator}/{R_c.denominator}")
    print(f"R(s) = 1 / tanh(ln(2)) = {R_s.numerator}/{R_s.denominator}")
    
    print("\n--- Final Equation ---")
    # To show the equation clearly, we print the components
    print(f"Value = {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator})) ^ ({m_0})")
    print(f"Value = {n_0} * ({ratio_R.numerator}/{ratio_R.denominator}) ^ ({m_0})")
    print(f"Value = {n_0} * ({ratio_R.denominator}/{ratio_R.numerator})")
    print(f"Value = {final_value.numerator}/{final_value.denominator}")
    
    print("\n--- Final Answer ---")
    print(f"{final_value.numerator}/{final_value.denominator}")


solve_fin_problem()
<<<252/17>>>