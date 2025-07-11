import math

def solve_triangle_median_range():
    """
    Calculates the range for the median m for angle A to be acute.

    The derivation shows that m is related to an intermediate variable y (where y=k^2)
    by the formula m^2 = (720 / (y - 25))^2 + 144.
    The condition for angle A to be acute constrains y to the union of two intervals:
    (0, 25) and (84.5, 169).

    This function calculates the values of m at the boundaries of these intervals.
    """

    print("Step 1: Calculate the boundaries for the first interval of m.")
    # This interval for m corresponds to y in (84.5, 169).
    # The lower boundary for m occurs as y approaches 169.
    y_lim1 = 169.0
    f_at_y1 = 720.0 / (y_lim1 - 25.0)  # f = 720 / 144 = 5
    m_sq_at_y1 = f_at_y1**2 + 144.0
    m_at_y1 = math.sqrt(m_sq_at_y1) # m = sqrt(25 + 144) = 13

    # The upper boundary for m occurs as y approaches 84.5.
    y_lim2 = 84.5
    f_at_y2 = 720.0 / (y_lim2 - 25.0)  # f = 720 / 59.5 = 1440/119
    m_sq_at_y2 = f_at_y2**2 + 144.0
    m_at_y2 = math.sqrt(m_sq_at_y2) # m = sqrt((1440/119)^2 + 144) = 2028/119
    m_at_y2_frac_num = 2028
    m_at_y2_frac_den = 119
    
    print("Step 2: Calculate the boundaries for the second interval of m.")
    # This interval for m corresponds to y in (0, 25).
    # The lower boundary for m occurs as y approaches 0.
    y_lim3 = 0.0
    f_at_y3 = 720.0 / (y_lim3 - 25.0) # f = 720 / -25 = -28.8
    m_sq_at_y3 = f_at_y3**2 + 144.0
    m_at_y3 = math.sqrt(m_sq_at_y3) # m = sqrt((-28.8)^2 + 144) = 31.2
    m_at_y3_frac_num = 156
    m_at_y3_frac_den = 5
    
    # The upper boundary for m occurs as y approaches 25, which tends to infinity.
    
    print("\nFinal Result:")
    print("The analysis reveals two disjoint intervals for the value of m.")
    print("The final range for m is given by the following inequalities:")
    print(f"{int(m_at_y1)} < m < {m_at_y2_frac_num}/{m_at_y2_frac_den} or m > {m_at_y3_frac_num}/{m_at_y3_frac_den}")

solve_triangle_median_range()