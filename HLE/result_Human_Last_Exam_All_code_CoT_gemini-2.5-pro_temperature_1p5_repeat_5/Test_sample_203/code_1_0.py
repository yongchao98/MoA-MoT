import math

def solve_median_range():
    """
    Calculates the range of the median m for which angle A is acute.
    """
    # Given values
    ad = 12  # Height
    ae = 13  # Angle bisector

    # From right triangle ADE, we find DE
    # de^2 = ae^2 - ad^2
    de_sq = ae**2 - ad**2
    de = math.sqrt(de_sq)

    # From the angle bisector theorem and geometric properties, we derive the relation:
    # (ad^2 - de^2) * 2 * DF = 2 * de * (ad^2 + BD*DC)
    # (144 - 25) * DF = 5 * (144 + BD*DC)
    # 119 * DF = 5 * (144 + BD*DC)
    # where DF is the distance from D to the median's base F.
    
    # The condition for angle A to be acute is BD*DC < AD^2 = 144.
    # For a valid triangle, BD > 0 and DC > 0, so BD*DC > 0.
    
    # This gives us a range for DF:
    # Lower bound for DF (when BD*DC -> 0):
    df_lower = (5 * (ad**2 + 0)) / (ad**2 - de_sq)
    # Upper bound for DF (when BD*DC -> 144):
    df_upper = (5 * (ad**2 + ad**2)) / (ad**2 - de_sq)
    
    # Now we find the range for m using the relation m^2 = DF^2 + AD^2
    m_sq_lower = df_lower**2 + ad**2
    m_sq_upper = df_upper**2 + ad**2
    
    m_lower = math.sqrt(m_sq_lower)
    m_upper = math.sqrt(m_sq_upper)
    
    # For the exact fraction form of the answer
    # Lower bound of m:
    # m_lower^2 = (720/119)^2 + 144 = (12^2 * 60^2 / 119^2) + 12^2
    # m_lower^2 = (12^2/119^2) * (60^2 + 119^2)
    # 60^2 + 119^2 = 3600 + 14161 = 17761
    num1 = 12
    num2 = 119**2 + 60**2
    den = 119
    
    # Upper bound of m:
    # m_upper^2 = (1440/119)^2 + 144 = (12^2 * 120^2 / 119^2) + 12^2
    # m_upper^2 = (12^2/119^2) * (120^2 + 119^2)
    # 120^2 + 119^2 = 14400 + 14161 = 28561 = 169^2
    num3 = 12 * 169
    
    lower_bound_expr = f"({num1} * sqrt({num2})) / {den}"
    upper_bound_expr = f"{num3} / {den}"

    print("The range of values for m is between the lower and upper bounds.")
    print(f"Lower bound: {lower_bound_expr}")
    print(f"Upper bound: {upper_bound_expr}")
    print("\nNumerically, the range is approximately:")
    print(f"({m_lower:.4f}, {m_upper:.4f})")

solve_median_range()
<<<The range of values for m is ((12 * sqrt(17761)) / 119, 2028 / 119)>>>