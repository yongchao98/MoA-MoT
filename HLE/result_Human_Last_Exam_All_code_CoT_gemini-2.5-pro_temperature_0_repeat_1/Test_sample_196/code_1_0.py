import sympy
from sympy import pi, Symbol, Rational, integrate

def solve_volume_problem():
    """
    Calculates the volume of the space enclosed by a cone and an ellipsoid.
    """
    # --- Step 1: Define parameters and find the intersection ---
    # Ellipsoid S2: x^2/3 + y^2/4 + z^2/3 = 1
    # Cone S1: Vertex at (0, 4, 0), tangent to S2.
    # The plane of tangency is found to be y=1.
    # Substituting y=1 into the ellipsoid equation gives the intersection circle:
    # x^2/3 + 1/4 + z^2/3 = 1  =>  x^2 + z^2 = 9/4
    R = Rational(3, 2)  # Radius of the intersection circle
    y_intersect = 1
    
    print("The cone and ellipsoid are tangent along a circle of radius R = 3/2 in the plane y = 1.")
    print("-" * 50)

    # --- Step 2: Calculate the volume of the cone part (V_cone) ---
    # This is a cone with vertex at y=4 and base at y=1.
    h_cone = 4 - y_intersect
    
    # Using the formula V = (1/3) * pi * R^2 * h
    V_cone = Rational(1, 3) * pi * (R**2) * h_cone
    
    print("Step 1: Calculate the volume of the cone from its vertex (y=4) to the tangency plane (y=1).")
    print(f"The cone's height is h = 4 - {y_intersect} = {h_cone}.")
    print(f"The cone's base radius is R = {R}.")
    print(f"V_cone = (1/3) * pi * R^2 * h")
    print(f"V_cone = (1/3) * pi * ({R})^2 * {h_cone}")
    print(f"V_cone = (1/3) * pi * {R**2} * {h_cone} = {V_cone}")
    print("-" * 50)

    # --- Step 3: Calculate the volume of the ellipsoid cap (V_cap) ---
    # This is the volume of the ellipsoid for y from 1 to its top at y=2.
    y = Symbol('y')
    a_sq = 3
    b_sq = 4
    
    # The cross-sectional area A(y) = pi * r(y)^2
    # From the ellipsoid equation, r(y)^2 = a^2 * (1 - y^2/b^2)
    r_sq_ellipsoid = a_sq * (1 - y**2 / b_sq)
    integrand = pi * r_sq_ellipsoid
    
    # Integrate A(y) from y=1 to y=2
    V_cap = integrate(integrand, (y, 1, 2))
    
    print("Step 2: Calculate the volume of the ellipsoid cap above the tangency plane (y=1).")
    print("The volume is found by integrating the cross-sectional area A(y) from y=1 to y=2.")
    print(f"A(y) = pi * r(y)^2 = pi * {a_sq} * (1 - y^2/{b_sq})")
    print(f"V_cap = integral from 1 to 2 of ({integrand}) dy")
    
    # Show the steps of integration
    antiderivative = integrate(integrand, y)
    val_at_2 = antiderivative.subs(y, 2)
    val_at_1 = antiderivative.subs(y, 1)
    
    print(f"V_cap = [ {antiderivative} ] from 1 to 2")
    print(f"V_cap = ( {val_at_2} ) - ( {val_at_1} )")
    print(f"V_cap = {V_cap}")
    print("-" * 50)

    # --- Step 4: Calculate the final enclosed volume ---
    V_total = V_cone - V_cap
    
    print("Step 3: Calculate the total enclosed volume.")
    print("The total volume is the difference between the cone's volume and the ellipsoid cap's volume.")
    print(f"V_total = V_cone - V_cap")
    print(f"V_total = {V_cone} - {V_cap} = {V_total}")
    
    # Print the final answer in the required format
    print("\nFinal Answer:")
    print(f"The final volume is {V_total.evalf()}.")
    
    return V_total.evalf()

if __name__ == '__main__':
    final_answer = solve_volume_problem()
    print(f"<<<{final_answer}>>>")
