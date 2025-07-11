import sympy

def solve_volume():
    """
    Calculates the volume of the space enclosed by a cone and an ellipsoid.
    """
    # Define symbols for symbolic math
    y = sympy.Symbol('y')
    pi = sympy.pi

    # --- Step 1: Calculate the volume of the cone section (V_cone) ---
    # The cone has its vertex at (0, 4, 0) and its base on the plane y=1.
    # The height of this cone is h = 4 - 1 = 3.
    # The radius of the base (from x^2 + z^2 = 9/4 at y=1) is r = 3/2.
    h_cone = 3
    r_base = sympy.Rational(3, 2)

    # Volume formula for a cone: (1/3) * pi * r^2 * h
    V_cone = (sympy.Rational(1, 3)) * pi * r_base**2 * h_cone
    print(f"Volume of the cone (V_cone) from y=1 to y=4 is calculated as (1/3) * pi * (3/2)^2 * 3.")
    print(f"V_cone = {V_cone}\n")

    # --- Step 2: Calculate the volume of the ellipsoid cap (V_cap) ---
    # We integrate the cross-sectional area of the ellipsoid from y=1 to y=2.
    # The ellipsoid is x^2/3 + y^2/4 + z^2/3 = 1.
    # A cross-section at height y is a circle x^2 + z^2 = 3 * (1 - y^2/4).
    # The area of this circle is A(y) = pi * (radius)^2 = pi * 3 * (1 - y^2/4).
    area_cross_section = pi * 3 * (1 - y**2/4)

    # Integrate the area from y=1 to y=2.
    V_cap = sympy.integrate(area_cross_section, (y, 1, 2))
    print("Volume of the ellipsoid cap (V_cap) from y=1 to y=2 is the integral of the cross-sectional area A(y) = pi * 3 * (1 - y^2/4).")
    print(f"V_cap = Integral(A(y) dy) from 1 to 2 = {V_cap}\n")

    # --- Step 3: Calculate the final enclosed volume ---
    # The enclosed volume is the difference between the cone's volume and the ellipsoid cap's volume.
    V_final = V_cone - V_cap
    
    print("The final enclosed volume is the difference between the volume of the cone and the volume of the ellipsoid cap.")
    # The problem asks to output the numbers in the final equation.
    print(f"Final Volume = V_cone - V_cap")
    print(f"Final Volume = {V_cone} - {V_cap} = {V_final}")


solve_volume()