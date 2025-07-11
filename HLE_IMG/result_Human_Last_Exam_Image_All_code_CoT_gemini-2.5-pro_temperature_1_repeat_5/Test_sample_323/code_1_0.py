import math
from fractions import Fraction

def solve_heat_transfer_problem():
    """
    Solves the heat transfer problem by calculating the required values step-by-step.
    """
    # Step 1: Determine n_0 and m_0
    # From the problem description and analysis of heat transfer principles (Q proportional to sqrt(k)),
    # and the Q values from the plots:
    # Q1=11.6, Q3=7.6, Q5=8.4, Q7=12.7, Q9=7.19
    # The lowest Q corresponds to Carbon Steel.
    n_0 = 9

    # Based on ratio analysis (as detailed in the thinking steps), the geometry of the fin in plot n_0=9
    # is determined to be circular.
    # m_0 = 1 for circular, m_0 = -1 for square.
    m_0 = 1

    print(f"Step 1: Determine n_0 and m_0")
    print(f"The carbon steel fin corresponds to the plot with the lowest heat transfer, which is plot {n_0}.")
    print(f"The geometry of fin {n_0} is determined to be circular.")
    print(f"Therefore, n_0 = {n_0} and m_0 = {m_0}\n")

    # Step 2: Calculate R(c) and R(s)

    # For the circular fin:
    # Given: Lh/k = ln(13), 4L/d = ln(13)
    # m^2 * L^2 = (4h/(kd)) * L^2 = (4hL/k) * (L/d) = 4 * ln(13) * (ln(13)/4) = (ln(13))^2
    # So, mL = ln(13)
    # And h/(mk) = (hL/k) / (mL) = ln(13) / ln(13) = 1
    M_c = math.log(13)
    N_c = 1

    # tanh(ln(x)) = (x - 1/x) / (x + 1/x) = (x^2 - 1) / (x^2 + 1)
    # coth(ln(x)) = (x^2 + 1) / (x^2 - 1)
    tanh_Mc = (13**2 - 1) / (13**2 + 1) # 168 / 170 = 84 / 85
    coth_Mc = (13**2 + 1) / (13**2 - 1) # 170 / 168 = 85 / 84

    R_c_num = 1 + N_c * coth_Mc
    R_c_den = 1 + N_c * tanh_Mc
    
    # Using Fractions for precision
    R_c = Fraction(1 + Fraction(85, 84)) / Fraction(1 + Fraction(84, 85))
    
    print(f"Step 2: Calculate R(c) and R(s)")
    print(f"For the circular fin, R(c) is calculated.")
    print(f"R(c) = (1 + coth(ln(13))) / (1 + tanh(ln(13))) = {R_c.numerator}/{R_c.denominator}")

    # For the square fin:
    # Given: Lh/k = ln(2), 4L/w = ln(2)
    # mL = ln(2), h/mk = 1
    M_s = math.log(2)
    N_s = 1

    tanh_Ms = (2**2 - 1) / (2**2 + 1) # 3 / 5
    coth_Ms = (2**2 + 1) / (2**2 - 1) # 5 / 3

    # Using Fractions for precision
    R_s = Fraction(1 + Fraction(5, 3)) / Fraction(1 + Fraction(3, 5))

    print(f"For the square fin, R(s) is calculated.")
    print(f"R(s) = (1 + coth(ln(2))) / (1 + tanh(ln(2))) = {R_s.numerator}/{R_s.denominator}\n")
    
    # Step 3: Final Calculation
    # expression = n_0 * (R(c)/R(s))^m_0
    
    ratio_R = R_c / R_s
    final_result = n_0 * (ratio_R ** m_0)
    
    print("Step 3: Final Calculation")
    print("The final expression is n_0 * (R(c)/R(s))^m_0")
    print(f"Substituting the values:")
    # "Remember in the final code you still need to output each number in the final equation!"
    print(f"Equation: {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^{m_0}")
    print(f"Result = {n_0} * ({ratio_R.numerator}/{ratio_R.denominator})^{m_0}")
    print(f"Result = {final_result.numerator}/{final_result.denominator}")

solve_heat_transfer_problem()