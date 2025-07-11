import math

def solve_triangle_median_range():
    """
    Calculates the range for the median m based on the derived inequality.
    """
    # The derived inequality for m^2 is:
    # 144 + (720/119)^2 < m^2 < 144 + (1440/119)^2

    # Calculate the lower bound for m
    lower_bound_sq = 144 + (720 / 119)**2
    lower_bound_m = math.sqrt(lower_bound_sq)

    # Calculate the upper bound for m
    upper_bound_sq = 144 + (1440 / 119)**2
    upper_bound_m = math.sqrt(upper_bound_sq)
    
    # We can simplify the upper bound expression algebraically:
    # m^2 < 144 + (1440/119)^2 = 144 * (1 + 14400 / (144 * 119^2)) = 144 * (1 + 100/119^2) = 144 * (119^2 + 100) / 119^2 is wrong
    # m^2 < 144 + (1440/119)^2 = 144 + 2073600 / 14161 = (144*14161 + 2073600)/14161 = 4112784/14161
    # Note that 4112784 = (2028)^2, so m^2 < (2028/119)^2
    upper_bound_m_fraction = "2028/119"

    # The lower bound doesn't simplify to a rational number.
    # lower_bound_sq = 144 + 518400/14161 = 2557584/14161
    lower_bound_m_fraction = "sqrt(2557584)/119"
    
    print("The range of values for m such that angle A is acute is given by the inequality:")
    print(f"   {lower_bound_m_fraction} < m < {upper_bound_m_fraction}")
    print("\nWhich corresponds to the numerical range:")
    print(f"   {lower_bound_m:.4f} < m < {upper_bound_m:.4f}")

solve_triangle_median_range()