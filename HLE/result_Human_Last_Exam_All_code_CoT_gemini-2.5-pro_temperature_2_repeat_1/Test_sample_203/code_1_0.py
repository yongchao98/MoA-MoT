import math

def solve_for_m_range():
    """
    This function calculates the range of the median m for angle A to be acute,
    based on the derivation explained above.
    """
    # From the derivation, we have the inequality for m^2:
    # m^2 < 144 + (1440 / 119)^2
    
    # Let's calculate the value of the upper bound.
    # m_upper^2 = 144 + (1440/119)^2
    # m_upper = sqrt(12^2 + (12*120/119)^2)
    # m_upper = 12 * sqrt(1 + (120/119)^2)
    # m_upper = 12 * sqrt( (119^2 + 120^2) / 119^2 )
    # m_upper = (12/119) * sqrt(14161 + 14400)
    # m_upper = (12/119) * sqrt(28561)
    
    # We find that 28561 is a perfect square.
    sqrt_28561 = int(math.sqrt(28561)) # This is 169
    
    # So the upper bound is m < (12 * 169) / 119
    m_upper_bound_numerator = 12 * sqrt_28561
    m_upper_bound_denominator = 119

    # The lower bound for m is 13.
    lower_bound = 13

    # Print the final result in the required format.
    # The final equation is the inequality that defines the range of m.
    print("The range of values for m for which angle A is acute is given by the inequality:")
    print(f"{lower_bound} < m < {m_upper_bound_numerator}/{m_upper_bound_denominator}")
    
solve_for_m_range()