import sympy
from sympy import pi, pprint

def solve_volume_problem():
    """
    Calculates the volume enclosed by a cone S1 and an ellipsoid S2.
    The code explains each step of the calculation.
    """
    # Define the symbolic variable for integration
    y = sympy.Symbol('y')

    # Step 1: Define the surfaces and find their intersection.
    # The ellipsoid S2 is given by x^2/3 + y^2/4 + z^2/3 = 1.
    # The cone S1 has its vertex at (0, 4, 0) and is tangent to S2.
    # Using the polar plane formula, the plane of tangency is found to be y = 1.
    # The volume is bounded by this plane and the top of the ellipsoid (at y=2).
    lower_limit = 1
    upper_limit = 2

    print("Step 1: Finding the boundaries for the integration.")
    print("The cone S1 is tangent to the ellipsoid S2.")
    print("The plane of tangency is calculated to be y = 1.")
    print("The ellipsoid S2 extends up to y = 2.")
    print(f"Thus, we will integrate along the y-axis from y = {lower_limit} to y = {upper_limit}.\n")

    # Step 2: Determine the cross-sectional areas.
    # For any y, the cross-sections are circles. We need their radii squared (r^2).

    # For the ellipsoid S2: x^2/3 + z^2/3 = 1 - y^2/4  =>  x^2 + z^2 = 3*(1 - y^2/4)
    r_e_sq = 3 * (1 - y**2 / 4)

    # For the cone S1: Its radius is 0 at its vertex y=4.
    # At the tangency plane y=1, its radius matches the ellipsoid's.
    # r^2 at y=1 is 3*(1 - 1^2/4) = 9/4, so the radius R = 3/2.
    # The height of this cone section is h = 4 - 1 = 3.
    # By similar triangles, the radius r_c at height y is r_c/(4-y) = R/h = (3/2)/3 = 1/2.
    # So, r_c = (4-y)/2.
    r_c_sq = ((4 - y) / 2)**2

    print("Step 2: Defining the cross-sectional areas.")
    print(f"Squared radius of the cone's cross-section: r_c^2 = {r_c_sq}")
    print(f"Squared radius of the ellipsoid's cross-section: r_e^2 = {r_e_sq}")

    # The area of the annular region at height y is A(y) = pi * (r_c^2 - r_e^2).
    integrand = pi * (r_c_sq - r_e_sq)
    simplified_integrand = sympy.simplify(integrand)
    print(f"The integrand A(y) = pi * (r_c^2 - r_e^2) simplifies to: {simplified_integrand}\n")

    # Step 3: Set up and evaluate the volume integral.
    volume_integral = sympy.Integral(simplified_integrand, (y, lower_limit, upper_limit))
    print("Step 3: Setting up the volume integral.")
    print("Volume V is the integral of A(y) from y=1 to y=2:")
    pprint(volume_integral)
    print()

    # To show the numbers in the final equation, we find the antiderivative first.
    antiderivative = sympy.integrate(simplified_integrand, y)
    print(f"The antiderivative of {simplified_integrand} is: {antiderivative}\n")

    print("Step 4: Evaluating the definite integral.")
    print("V = [antiderivative evaluated at y=2] - [antiderivative evaluated at y=1]")
    print(f"V = [{antiderivative}] from y=1 to y=2")
    # Show the numerical substitution as requested
    print(f"V = pi * [ (({upper_limit})**3/3 - ({upper_limit})**2 + ({upper_limit})) - (({lower_limit})**3/3 - ({lower_limit})**2 + ({lower_limit})) ]")
    term_upper = f"{upper_limit**3}/{3} - {upper_limit**2} + {upper_limit}"
    term_lower = f"{lower_limit**3}/{3} - {lower_limit**2} + {lower_limit}"
    print(f"V = pi * [ ({term_upper}) - ({term_lower}) ]")
    val_upper = (upper_limit**3/3 - upper_limit**2 + upper_limit)
    val_lower = (lower_limit**3/3 - lower_limit**2 + lower_limit)
    print(f"V = pi * [ ({val_upper:.2f}) - ({val_lower:.2f}) ]")
    print(f"V = pi * [ (8/3 - 2) - (1/3) ]")
    print(f"V = pi * [ (2/3) - (1/3) ]")
    print(f"V = pi * (1/3)\n")
    
    # Calculate the final volume
    volume = sympy.integrate(simplified_integrand, (y, lower_limit, upper_limit))
    
    print("Final Result:")
    print(f"The volume of the space enclosed by the two surfaces is {volume}.")

if __name__ == "__main__":
    solve_volume_problem()