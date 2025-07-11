import math

def solve_triangle_median_range():
    """
    Calculates the range of the median m for which angle A is acute.
    
    The derivation uses coordinate geometry, leading to the following key relationships:
    - Let u = y - x and v = xy, where B=(-x,0) and C=(y,0).
    - Condition from angle bisector: 119u - 10v = 1440
    - Median length: m^2 = u^2/4 + 144
    - Condition for angle A to be acute: 0 < v < 144
    """
    
    # 1. Find the range for u from the range of v (0 < v < 144)
    # u = (1440 + 10v) / 119
    v_lower = 0
    v_upper = 144
    
    u_lower = (1440 + 10 * v_lower) / 119
    u_upper = (1440 + 10 * v_upper) / 119
    
    # 2. Find the range for m from the range of u
    # m = sqrt(u^2/4 + 144)
    
    # Calculate lower bound for m
    m_sq_lower = (u_lower / 2)**2 + 144
    m_lower = math.sqrt(m_sq_lower)
    
    # Calculate upper bound for m
    m_sq_upper = (u_upper / 2)**2 + 144
    m_upper = math.sqrt(m_sq_upper)
    
    print("The range for the median m is (m_lower, m_upper), where:")
    # Using fractions for exact representation before taking the root
    from fractions import Fraction
    # For upper bound: m_upper = sqrt( ( (2880/119)/2 )^2 + 144 ) = sqrt( (1440/119)^2 + 144)
    # m_upper^2 = (1440^2 + 144 * 119^2) / 119^2 = (144*10^2 + 144*119^2)/119^2... no
    # m_upper^2 = 144 * ((120/119)^2 + 1) = 144 * (120^2 + 119^2) / 119^2 = 144 * 28561 / 119^2
    # sqrt(28561)=169. So m_upper = 12*169/119 = 2028/119
    m_upper_exact_num = 2028
    m_upper_exact_den = 119

    # For lower bound: m_lower = sqrt( ( (1440/119)/2 )^2 + 144 ) = sqrt( (720/119)^2 + 144 )
    # m_lower^2 = 144 * ((60/119)^2+1) = 144 * (60^2 + 119^2)/119^2 = 144 * 17761 / 119^2
    # so m_lower = 12 * sqrt(17761) / 119
    m_lower_exact_term1 = 12
    m_lower_exact_term2 = 17761
    m_lower_exact_den = 119
    
    print(f"The lower bound m_lower = ({m_lower_exact_term1} * sqrt({m_lower_exact_term2})) / {m_lower_exact_den}")
    print(f"Which is approximately: {m_lower}")
    
    print(f"\nThe upper bound m_upper = {m_upper_exact_num} / {m_upper_exact_den}")
    print(f"Which is approximately: {m_upper}")

solve_triangle_median_range()