import math

def solve_triangle_problem():
    """
    Calculates the range of the median m for angle A to be acute.
    """
    # Lower bound calculation
    # m^2 > 144 + (720/119)^2
    m_sq_lower_num = 144 * 119**2 + 720**2
    m_sq_lower_den = 119**2
    
    # The numerator is 2557584 = 144 * 17761
    # So the lower bound for m is sqrt(144 * 17761) / 119
    lower_bound_val_1 = 12
    lower_bound_val_2 = 17761
    lower_bound_den = 119
    
    # Upper bound calculation
    # m^2 < 144 + (1440/119)^2
    m_sq_upper_num = 144 * 119**2 + 1440**2
    m_sq_upper_den = 119**2
    
    # The numerator is 4112784 = 2028^2
    # So the upper bound for m is 2028 / 119
    upper_bound_num = int(math.sqrt(m_sq_upper_num))
    upper_bound_den = 119

    # Calculate numerical values for printing
    lower_bound_numerical = (lower_bound_val_1 * math.sqrt(lower_bound_val_2)) / lower_bound_den
    upper_bound_numerical = upper_bound_num / upper_bound_den

    print("The range of values for m is given by the inequality:")
    print(f"({lower_bound_val_1} * sqrt({lower_bound_val_2})) / {lower_bound_den} < m < {upper_bound_num} / {upper_bound_den}")
    print("\nNumerically, this is approximately:")
    print(f"{lower_bound_numerical:.4f} < m < {upper_bound_numerical:.4f}")

solve_triangle_problem()