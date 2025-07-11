import math

def solve_geometry_problem():
    """
    Calculates the range of values for the median m for which angle A is acute.
    """
    # Given values
    AD = 12
    AE = 13

    # From Pythagorean theorem in triangle ADE, DE = 5.
    DE = math.sqrt(AE**2 - AD**2)

    # The problem can be solved using coordinate geometry, which leads to two ranges for m.
    # The boundaries of these ranges are determined by specific configurations of the triangle.

    # Boundary 1: A degenerate triangle where B, C, and E coincide.
    # In this case, the median AF is the same as the angle bisector AE.
    m_lower_bound1 = AE
    
    # Boundary 2: The case where angle A is a right angle (90 degrees).
    # The analysis leads to k = 84.5 = 169/2.
    # m^2 = (-720 / (25 - 169/2))^2 + 144
    # m = 12 * 169 / 119
    m_upper_bound1_num = 12 * 169
    m_upper_bound1_den = 119
    m_upper_bound1_val = m_upper_bound1_num / m_upper_bound1_den

    # Boundary 3: Another degenerate triangle case.
    # The analysis leads to k = 0.
    # m^2 = (-720 / (25 - 0))^2 + 144
    m_lower_bound2_val = math.sqrt((-720 / 25)**2 + 144)
    # This value is 31.2 or 156/5
    m_lower_bound2_num = 156
    m_lower_bound2_den = 5

    # The final range for m is the union of two intervals.
    print("For angle A to be acute, the range of values for the median m is the union of two intervals.")
    print(f"The first interval is: {m_lower_bound1} < m < {m_upper_bound1_num}/{m_upper_bound1_den} (approx. {m_upper_bound1_val:.3f})")
    print(f"The second interval is: m > {m_lower_bound2_num}/{m_lower_bound2_den} (which is {m_lower_bound2_val})")
    print("\nFinal Answer:")
    print(f"13 < m < 2028/119  or  m > 156/5")

solve_geometry_problem()