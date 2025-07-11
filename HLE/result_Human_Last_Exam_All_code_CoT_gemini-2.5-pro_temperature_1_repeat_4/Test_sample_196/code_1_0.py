import sympy as sp

def find_enclosed_volume():
    """
    This script calculates the volume of the space enclosed by a cone S1 and an ellipsoid S2.

    S1: Cone with vertex at (0, 4, 0), tangent to S2.
    S2: Ellipsoid x^2/3 + y^2/4 + z^2/3 = 1.

    The method involves finding the plane of tangency and then integrating the difference
    in the cross-sectional areas of the two surfaces.
    """
    # Define symbols for our calculations
    x, y, z, x0, y0, z0 = sp.symbols('x y z x0 y0 z0')

    # Ellipsoid S2 parameters
    a_sq, b_sq, c_sq = 3, 4, 3
    
    # Cone S1 vertex
    V_vertex = (0, 4, 0)

    # --- Step 1: Find the plane of tangency ---
    print("--- Step 1: Find the plane of tangency ---")
    print(f"The ellipsoid S2 is x^2/{a_sq} + y^2/{b_sq} + z^2/{c_sq} = 1.")
    print(f"The cone S1 has its vertex at {V_vertex}.")
    
    # The tangent plane to the ellipsoid at a point (x0, y0, z0) is (x*x0)/a^2 + (y*y0)/b^2 + (z*z0)/c^2 = 1
    # This plane must contain the cone's vertex. We substitute the vertex coordinates into the plane equation.
    # (0*x0)/3 + (4*y0)/4 + (0*z0)/3 = 1  =>  y0 = 1
    y_tangent = 1
    print(f"By substituting the cone's vertex into the general tangent plane equation, we find the plane of tangency is y = {y_tangent}.")
    print("The volume is enclosed between this plane and the top of the ellipsoid.\n")

    # --- Step 2: Define the cross-sectional areas A(y) ---
    print("--- Step 2: Define the cross-sectional areas A(y) ---")
    # For the ellipsoid S2, the cross-section is a circle with radius r_e^2 = 3*(1 - y^2/4)
    r_sq_ellipsoid = a_sq * (1 - y**2/b_sq)
    Area_ellipsoid = sp.pi * r_sq_ellipsoid
    print(f"Ellipsoid cross-sectional area A_ellipsoid(y) = pi * {sp.simplify(r_sq_ellipsoid)}")
    
    # For the cone S1, its radius is r=0 at the vertex y=4.
    # At the tangency plane y=1, its radius matches the ellipsoid's.
    # r_c(y=1)^2 = 3*(1 - 1^2/4) = 9/4, so r_c(y=1) = 3/2.
    # The cone's radius r_c as a function of y is linear: r_c(y) = (4-y)/2.
    r_sq_cone = ((V_vertex[1] - y)/2)**2
    Area_cone = sp.pi * r_sq_cone
    print(f"Cone cross-sectional area A_cone(y) = pi * {sp.simplify(r_sq_cone)}\n")

    # --- Step 3: Set up and solve the volume integral ---
    print("--- Step 3: Calculate the volume by integration ---")
    y_top = sp.sqrt(b_sq) # Top of ellipsoid is at y=2
    integrand = Area_cone - Area_ellipsoid
    
    # The integrand simplifies to pi * (y-1)^2
    integrand_simplified = sp.simplify(integrand)
    
    print(f"The volume V is the integral of (A_cone(y) - A_ellipsoid(y)) from y = {y_tangent} to y = {y_top}.")
    print(f"V = Integral from {y_tangent} to {y_top} of ({integrand_simplified}) dy")
    
    # Show the integration steps
    antiderivative = sp.integrate(integrand_simplified/sp.pi, y)
    print(f"\nFirst, we find the antiderivative of the simplified integrand (y-1)^2:")
    print(f"Integral[(y-1)^2] dy = {antiderivative}")
    
    print("\nNext, we evaluate this antiderivative at the limits of integration, y=2 and y=1:")
    print(f"V = pi * [ {antiderivative} ] from 1 to 2")
    
    upper_bound_val = antiderivative.subs(y, 2)
    lower_bound_val = antiderivative.subs(y, 1)

    print(f"V = pi * ( (({2}-1)^3/{3}) - (({1}-1)^3/{3}) )")
    print(f"V = pi * ( ({1**3}/{3}) - ({0**3}/{3}) )")
    print(f"V = pi * ( {sp.Rational(1,3)} - {0} )")
    
    # Final Result
    final_volume = sp.integrate(integrand, (y, y_tangent, y_top))
    print(f"\nThe calculated volume is {final_volume}.")

if __name__ == '__main__':
    find_enclosed_volume()
    # The final answer is pi/3.
    print("\n<<<pi/3>>>")