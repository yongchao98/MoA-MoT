import sympy as sp

def solve_volume_problem():
    """
    Solves the volume problem step-by-step and prints the explanation.
    The volume of the space enclosed by the cone S1 and the ellipsoid S2 is calculated.
    """
    # Define symbolic variable
    y = sp.symbols('y')

    # Introduction
    print("This script calculates the volume of the space enclosed by a cone S1 and an ellipsoid S2.")
    print("S1: Cone with vertex at (0,4,0), tangent to S2.")
    print("S2: Ellipsoid given by x^2/3 + y^2/4 + z^2/3 = 1.")
    print("-" * 50)

    # Step 1: Find the plane of tangency
    print("Step 1: Finding the plane of tangency.")
    print("The tangent plane to the ellipsoid at (x0, y0, z0) must contain the cone's vertex (0, 4, 0).")
    print("This condition leads to y0 = 1. So, the tangency occurs in the plane y = 1.")
    y_tangency = 1
    y_top = 2
    print(f"The enclosed volume lies between y = {y_tangency} and the top of the ellipsoid at y = {y_top}.\n")

    # Step 2: Define the cross-sectional areas
    print("Step 2: Defining the cross-sectional areas A(y).")
    # Ellipsoid Area
    # from x^2/3 + z^2/3 = 1 - y^2/4  =>  x^2 + z^2 = r^2 = 3*(1 - y^2/4)
    A_ellipsoid = sp.pi * (sp.S(3)/4) * (4 - y**2)
    print(f"Area of ellipsoid cross-section: A_ellipsoid(y) = {A_ellipsoid}")

    # Cone Area
    # At y=1, r_cone^2 is the same as r_ellipsoid^2: 3/4 * (4 - 1^2) = 9/4
    # The cone's radius r(y) is linear from r=0 at y=4. So r(y) = k(4-y).
    # r(1) = 3/2 => k(3) = 3/2 => k = 1/2.
    # r(y) = (4-y)/2
    A_cone = sp.pi * ((4 - y)/2)**2
    print(f"Area of cone cross-section:      A_cone(y) = {sp.simplify(A_cone)}")
    print("-" * 50)

    # Step 3: Set up and solve the integral
    print("Step 3: Calculating the volume by integration.")
    print("V = Integral from 1 to 2 of [A_cone(y) - A_ellipsoid(y)] dy")

    integrand = A_cone - A_ellipsoid
    integrand_simplified = sp.simplify(integrand)
    print(f"\nThe integrand is: {integrand_simplified}")

    # Use the simplified form (y-1)^2 for clarity in the final calculation
    integrand_final_form = sp.pi * (y-1)**2
    antiderivative = sp.integrate(integrand_final_form, y)
    print(f"The indefinite integral is: {antiderivative}\n")

    # Step 4: Show the final calculation with numbers
    print("Step 4: The final equation and result.")
    # Values for substitution
    y_upper = 2
    y_lower = 1
    
    # Calculate terms for display
    upper_term = (y_upper - 1)**3
    lower_term = (y_lower - 1)**3
    
    print("Volume = [ pi/3 * (y - 1)^3 ] evaluated from y=1 to y=2.\n")
    print(f"Substituting the limits:")
    print(f"Volume = ( pi/3 * ({y_upper} - 1)^3 ) - ( pi/3 * ({y_lower} - 1)^3 )")
    print(f"Volume = ( pi/3 * ({upper_term})^3 ) - ( pi/3 * ({lower_term})^3 )")
    print(f"Volume = ( pi/3 * {upper_term} ) - ( pi/3 * {lower_term} )")
    print(f"Volume = {sp.pi/3} - 0")

    final_volume = sp.pi / 3
    print("\n" + "=" * 50)
    print(f"The final volume is {final_volume}")
    print(f"Numerically, the volume is approximately {final_volume.evalf():.4f}")
    print("=" * 50)

solve_volume_problem()