import math
from fractions import Fraction

def solve_heat_transfer_problem():
    """
    Solves the heat transfer problem by determining n_0, m_0, R(c), R(s)
    and then calculating the final expression.
    """
    # Part 1: Determine n_0 and m_0
    # Based on the analysis in the text, the carbon steel fin corresponds to plot 9,
    # which has the lowest heat transfer rate.
    n_0 = 9
    # The geometry of the fin in plot 9 is determined to be circular by comparing heat transfer
    # ratios and matching tip temperatures with the 3D models.
    # m_0 = 1 for circular geometry.
    m_0 = 1
    print(f"Step 1: Determine n_0 and m_0")
    print(f"The carbon steel fin is in plot n_0 = {n_0}")
    print(f"The geometry is circular, so m_0 = {m_0}")
    print("-" * 20)

    # Part 2: Calculate R(c) and R(s)
    print(f"Step 2: Calculate R(c) and R(s)")
    # For both cases, the given conditions lead to h/(m*k) = 1.
    # The ratio R simplifies to 1 / tanh(mL).

    # Circular case
    # mL = ln(13)
    # R(c) = 1 / tanh(ln(13))
    # tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85
    tanh_ln13_num = 13**2 - 1
    tanh_ln13_den = 13**2 + 1
    tanh_ln13 = Fraction(tanh_ln13_num, tanh_ln13_den)
    R_c = 1 / tanh_ln13
    print(f"For the circular fin, mL = ln(13) and h/mk = 1.")
    print(f"R(c) = 1 / tanh(ln(13)) = 1 / ({tanh_ln13.numerator}/{tanh_ln13.denominator}) = {R_c.numerator}/{R_c.denominator}")

    # Square case
    # mL = ln(2)
    # R(s) = 1 / tanh(ln(2))
    # tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5
    tanh_ln2_num = 2**2 - 1
    tanh_ln2_den = 2**2 + 1
    tanh_ln2 = Fraction(tanh_ln2_num, tanh_ln2_den) # This is (4-1)/(4+1) = 3/5
    R_s = 1 / tanh_ln2
    print(f"For the square fin, mL = ln(2) and h/mk = 1.")
    print(f"R(s) = 1 / tanh(ln(2)) = 1 / ({tanh_ln2.numerator}/{tanh_ln2.denominator}) = {R_s.numerator}/{R_s.denominator}")
    print("-" * 20)

    # Part 3: Final Calculation
    print(f"Step 3: Calculate the final expression")
    # Expression: n_0 * (R(c) / R(s)) ^ m_0
    ratio_R = R_c / R_s
    final_value = n_0 * (ratio_R) ** m_0

    print(f"The expression to calculate is: n_0 * (R(c) / R(s)) ^ m_0")
    print(f"Substituting the values: {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator})) ^ {m_0}")
    print(f"= {n_0} * ({ratio_R.numerator}/{ratio_R.denominator})")
    print(f"= ({n_0} * {ratio_R.numerator}) / {ratio_R.denominator}")
    print(f"= {final_value.numerator} / {final_value.denominator}")
    print("-" * 20)
    print(f"Final Answer: {final_value.numerator}/{final_value.denominator}")

solve_heat_transfer_problem()