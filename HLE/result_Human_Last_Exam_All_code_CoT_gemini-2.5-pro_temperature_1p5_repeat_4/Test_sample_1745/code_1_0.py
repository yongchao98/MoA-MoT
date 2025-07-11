import math

def display_potential_distribution():
    """
    This function prints the derived mathematical expression for the
    Electrical Double-Layer (EDL) potential distribution ψ(y).
    The expression is broken down into its symbolic components for clarity.
    """

    # --- Symbolic Representation of Variables ---
    # These strings represent the variables in the equation.
    potential_function = "ψ(y)"
    zeta_potential_1 = "z_1"
    slip_length = "β"
    debye_huckel_param = "k"
    channel_height = "H"
    vertical_position = "y"

    # --- Constructing the Final Equation String ---
    # The final equation is: ψ(y) = z_1 * (1 + β*k) * sinh(k*(H-y)) / sinh(k*H)
    
    # Term 1: Slip-dependent zeta potential at the bottom wall (z_a1)
    term1 = f"{zeta_potential_1} * (1 + {slip_length} * {debye_huckel_param})"
    
    # Term 2: Hyperbolic sine function component
    # Numerator: sinh(k*(H-y))
    # Denominator: sinh(k*H)
    term2_numerator = f"sinh({debye_huckel_param} * ({channel_height} - {vertical_position}))"
    term2_denominator = f"sinh({debye_huckel_param} * {channel_height})"
    
    # Full expression
    final_equation = f"{potential_function} = ({term1}) * ({term2_numerator} / {term2_denominator})"

    # --- Printing the Final Result ---
    print("The derived expression for the Electrical Double-Layer (EDL) potential distribution is:")
    print("-" * 70)
    print(final_equation)
    print("-" * 70)
    
    # As per the instructions, let's explicitly show the components that form the final equation
    print("\nWhere the equation components are:")
    print(f"  ψ(y): The EDL potential at vertical position y")
    print(f"  {zeta_potential_1}: The intrinsic zeta potential at the bottom surface (y=0)")
    print(f"  {slip_length}: The slip length")
    print(f"  {debye_huckel_param}: The Debye–Huckel parameter (inverse of EDL thickness)")
    print(f"  {channel_height}: The total height of the microchannel")
    print(f"  {vertical_position}: The vertical position within the channel (from 0 to H)")
    
display_potential_distribution()