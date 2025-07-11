import sympy as sp

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone.
    The cone has a base radius of 3 at z=0 and a height of 2 with its vertex at (0,0,2).
    """
    # Define the symbols for cylindrical coordinates
    r, theta, z = sp.symbols('r theta z')

    # Define the constants for the cone's geometry
    H = 2  # Height
    R = 3  # Radius

    # The function in cylindrical coordinates is f = z^2 * r^2
    # The volume element is dV = r dr d(theta) dz
    # The integrand is the product of the function and the volume element
    integrand = z**2 * r**3

    # Define the upper limit for the radius r at a given height z
    # The line from (r=R, z=0) to (r=0, z=H) is r(z) = R * (H - z) / H
    r_upper_limit = sp.S(R) * (H - z) / H

    # Perform the triple integration
    # Order of integration: r, then theta, then z
    integral_value = sp.integrate(
        integrand,
        (r, 0, r_upper_limit),
        (theta, 0, 2 * sp.pi),
        (z, 0, H)
    )

    # --- Output Section ---
    print("The integral to be solved is for the function f(x,y) = z^2 * (x^2 + y^2)")
    print("over a cone with height H=2 and base radius R=3.\n")
    print("In cylindrical coordinates, the setup for the final integral equation is:")
    
    # We explicitly print each number in the equation
    # The numbers are: z limits (0, 2), theta limits (0, 2*pi), r limits (0, 3/2*(2-z))
    # and the powers in the integrand (z^2, r^3).
    print(f"  Integral from z={0} to z={H} of (")
    print(f"    Integral from theta={0} to theta=2*pi of (")
    print(f"      Integral from r={0} to r=({R}/{H})*({H}-z) of ({z**2} * {r**3}) dr")
    print("    ) d(theta)")
    print("  ) dz")

    print("\nThe calculated exact value of the integral is:")
    print(integral_value)

solve_cone_integral()