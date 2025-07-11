import math

def solve_median_range():
    """
    Calculates the range of the median m based on the derived formulas.
    """
    # From the geometric analysis, we derived the range for a variable S = x1 - x2:
    # -2880/119 < S < -1440/119
    # And the median m is related to S by m^2 = S^2/4 + 144.

    # Calculate the bounds for S^2
    # Since S is in a negative interval, the function S^2 is decreasing over this interval.
    s_lower_bound_for_S_squared = (-1440 / 119)**2
    s_upper_bound_for_S_squared = (-2880 / 119)**2

    # Calculate the bounds for m^2
    m_squared_lower = s_lower_bound_for_S_squared / 4 + 144
    m_squared_upper = s_upper_bound_for_S_squared / 4 + 144

    # Calculate the bounds for m
    m_lower = math.sqrt(m_squared_lower)
    m_upper = math.sqrt(m_squared_upper)
    
    # We can simplify the upper bound calculation as it results in a known perfect square.
    # m_upper^2 = (1440/119)^2 + 144 = (144*10/119)^2 + 144 = 144 * (100/119^2 + 1)
    # = 144 * (100 + 119^2)/119^2 = 144 * (14400 + 14161)/14161 = 144 * 28561 / 14161
    # sqrt(28561) is 169.
    # So m_upper = 12 * 169 / 119 = 2028 / 119.
    
    m_upper_exact_val = 2028 / 119

    # The lower bound does not simplify to a rational number.
    # m_lower^2 = (720/119)^2 + 144 = 144 * ( (5*144/144)^2 / 119^2 + 1) -> No
    # m_lower^2 = (720^2 + 144 * 119^2) / 119^2 = (518400 + 144 * 14161) / 119^2
    # = (518400 + 2039184) / 14161 = 2557584 / 14161
    # = (144 * 17761) / 14161
    # m_lower = 12 * sqrt(17761) / 119
    
    # Printing the final equation with all its components.
    lower_bound_numerator = "12 * sqrt(17761)"
    denominator = "119"
    upper_bound_numerator = "2028"
    
    print(f"The range for m is given by the inequality:")
    print(f"{lower_bound_numerator} / {denominator} < m < {upper_bound_numerator} / {denominator}")
    print("\nWhich is approximately:")
    print(f"{m_lower:.3f} < m < {m_upper_exact_val:.3f}")

solve_median_range()