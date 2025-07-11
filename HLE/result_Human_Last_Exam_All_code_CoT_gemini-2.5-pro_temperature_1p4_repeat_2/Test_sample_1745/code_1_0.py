import math

def solve_and_print_expression():
    """
    This function derives and prints the final expression for the EDL potential distribution.
    The derivation is based on solving the linearized Poisson-Boltzmann equation
    with the given boundary conditions.
    """

    # Define the symbolic variables as strings to construct the final expression.
    # These represent the physical quantities in the problem.
    potential_distribution = "psi(y)"
    base_zeta_potential_1 = "z_1"
    slip_length = "beta"
    debye_huckel_parameter = "k"
    channel_height = "H"
    vertical_coordinate = "y"
    
    # The numbers present in the equation
    number_one = "1"
    number_two = "2"
    
    # The final expression is derived as:
    # psi(y) = z_1*(1 + beta*k) * (sinh(k*(H/2 - y)) / sinh(k*H))
    
    # We construct this expression using the string variables defined above.
    
    # The slip-dependent zeta potential part: z_1*(1 + beta*k)
    zeta_term = f"{base_zeta_potential_1}*({number_one} + {slip_length}*{debye_huckel_parameter})"
    
    # The hyperbolic sine ratio part: sinh(k*(H/2 - y)) / sinh(k*H)
    numerator_sinh = f"sinh({debye_huckel_parameter}*({channel_height}/{number_two} - {vertical_coordinate}))"
    denominator_sinh = f"sinh({debye_huckel_parameter}*{channel_height})"
    sinh_ratio = f"({numerator_sinh} / {denominator_sinh})"
    
    # Combine all parts into the final equation.
    final_equation = f"{potential_distribution} = {zeta_term} * {sinh_ratio}"
    
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    print(final_equation)
    
    print("\nWhere:")
    print(f"  psi(y): Electrical potential at vertical position y")
    print(f"  z_1: Base zeta potential at the bottom surface (y = -H/2)")
    print(f"  beta: Slip length")
    print(f"  k: Debye–Huckel parameter")
    print(f"  H: Height of the microchannel")
    print(f"  y: Vertical coordinate (-H/2 <= y <= H/2)")
    print(f"  sinh: Hyperbolic sine function")

# Execute the function to print the result.
solve_and_print_expression()
<<<psi(y) = z_1*(1 + beta*k) * (sinh(k*(H/2 - y)) / sinh(k*H))>>>