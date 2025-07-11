import math

def solve_triangle_median_range():
    """
    Calculates the range for the median m given the height AD=12,
    angle bisector AE=13, and that angle A is acute.
    """
    AD = 12

    # From the condition P = tan(B)tan(C) > 1, we found the range for DF.
    # The lower bound for DF corresponds to P -> infinity (A -> 0).
    # The upper bound for DF corresponds to P -> 1 (A -> 90 degrees).
    df_low = 720 / 119
    df_up = 1440 / 119

    # The median length m is related to AD and DF by m^2 = AD^2 + DF^2.
    # Since m = sqrt(AD^2 + DF^2) is an increasing function of DF for DF > 0,
    # we can find the range of m by using the bounds of DF.

    # Calculate the lower bound for m
    m_low_sq = AD**2 + df_low**2
    m_low = math.sqrt(m_low_sq)
    
    # Calculate the upper bound for m
    m_up_sq = AD**2 + df_up**2
    m_up = math.sqrt(m_up_sq)

    # The final equation for the range of m is m_low < m < m_up.
    # Let's find the exact fractional forms for the output.
    # m_up = sqrt(12^2 + (1440/119)^2) = 12/119 * sqrt(119^2 + 120^2)
    #      = (12/119) * sqrt(28561) = (12/119) * 169 = 2028/119
    m_up_num = 2028
    m_up_den = 119
    
    # m_low = sqrt(12^2 + (720/119)^2) = 12/119 * sqrt(119^2 + 60^2)
    #       = (12/119) * sqrt(17761)
    m_low_num_coef = 12
    m_low_sqrt_val = 17761
    m_low_den = 119

    print("The range of values for m is given by the inequality:")
    print(f"{m_low_num_coef}*sqrt({m_low_sqrt_val})/{m_low_den} < m < {m_up_num}/{m_up_den}")
    
solve_triangle_median_range()