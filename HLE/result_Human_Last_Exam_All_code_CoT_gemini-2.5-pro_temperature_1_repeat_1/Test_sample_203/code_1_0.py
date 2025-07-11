import math

def solve_triangle_median_range():
    """
    Calculates the range of the median m based on the given triangle properties.
    """
    # Given values
    h = 12  # AD, height
    l = 13  # AE, angle bisector length

    # Step 1: Calculate DE using Pythagorean theorem in triangle ADE
    # DE^2 = AE^2 - AD^2
    de_sq = l**2 - h**2
    de = math.sqrt(de_sq)

    # Step 2: Define key constants from the geometric derivation
    # The derivation leads to S_low = 10*h^2 / (h^2 - DE^2) and S_high = 2*S_low
    h_sq = h**2
    h_sq_minus_de_sq = h_sq - de_sq

    # Step 3: Calculate the upper bound for m (when angle A approaches 90 degrees)
    # The derivation shows m_high = h * l^2 / (h^2 - DE^2)
    m_high = h * (l**2) / h_sq_minus_de_sq

    # Step 4: Calculate the lower bound for m (when the triangle is right-angled at B or C)
    # The derivation shows m_low^2 = (5 * h^2 / (h^2 - DE^2))^2 + h^2
    s_low = (h_sq * 10) / h_sq_minus_de_sq
    m_low_sq = (s_low / 2)**2 + h_sq
    m_low = math.sqrt(m_low_sq)
    
    # Alternatively, using the simplified formula for m_low
    h_fourth = h_sq**2
    de_fourth = de_sq**2
    val_inside_sqrt = h_fourth - de_sq * h_sq + de_fourth
    m_low_alt = (h * math.sqrt(val_inside_sqrt)) / h_sq_minus_de_sq

    # The problem asks for the range of m
    # The lower bound is inclusive, the upper bound is exclusive.
    print(f"The height AD is {h}")
    print(f"The angle bisector AE is {l}")
    print(f"The median AF is m")
    print(f"For angle A to be acute, the range of values for m is:")
    print(f"{m_low_alt} <= m < {m_high}")
    # The final equation is the inequality showing the range of m.
    # We output each number in that final equation.
    print("\nFinal equation:")
    print(f"m >= {m_low_alt}")
    print(f"m < {m_high}")


solve_triangle_median_range()