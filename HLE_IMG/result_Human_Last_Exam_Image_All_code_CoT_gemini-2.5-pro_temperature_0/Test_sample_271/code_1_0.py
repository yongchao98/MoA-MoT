import math

def solve_equation():
    """
    This function defines the given parameters and prints the equation for v0
    with the numerical values substituted in.
    """
    # Given parameters from the problem description
    m1 = 1  # mass at the end of the rod (kg)
    m2 = 2  # mass of the projectile (kg)
    R = 3   # length of the rod (m)
    d = 1   # distance of impact from the pivot (m)
    g = 10  # acceleration due to gravity (m/s^2)

    # The derived equation for v0 is:
    # v0 = (2 / (m2 * d)) * sqrt(g * (m1*R + m2*d) * (m1*R^2 + m2*d^2))
    # We will now format this equation as a string with the given values.
    
    # Note: The Python ** operator is used for exponentiation (e.g., R**2 for R squared).
    
    equation_string = f"v0 = (2 / ({m2} * {d})) * sqrt({g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2))"

    print("The equation for the value that v0 must have is:")
    print(equation_string)

solve_equation()