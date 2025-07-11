import math

def display_final_equation():
    """
    This function prints the final derived expression for the Electrical
    Double-Layer (EDL) potential distribution psi(y).
    """

    # Define the components of the equation as strings for clear printing.
    # Note: Unicode characters are used for Greek letters for better readability.
    psi_y = "ψ(y)"
    z_1 = "z_1"
    beta = "β"
    k = "k"
    H = "H"
    y = "y"

    # The expression for the effective zeta potential at the bottom wall
    zeta_1_expr = f"{z_1} * (1 + {beta}*{k})"

    # The hyperbolic sine part of the equation
    # The numbers 2 are explicitly included as requested.
    numerator = f"sinh({k} * ({H}/2 - {y}))"
    denominator = f"sinh({k}*{H})"
    fraction_part = f"[{numerator} / {denominator}]"

    # Construct the full final equation string
    final_equation = f"{psi_y} = {zeta_1_expr} * {fraction_part}"

    print("The derived expression for the Electrical double-layer potential distribution is:")
    print(final_equation)
    print("\n--- Variable Definitions ---")
    print(f"ψ(y): Electrical potential at vertical position y.")
    print(f"{z_1}:   Zeta potential parameter for the bottom surface.")
    print(f"{beta}:   Slip length.")
    print(f"{k}:   Debye–Huckel parameter.")
    print(f"{H}:   Height of the microchannel.")
    print(f"{y}:   Vertical coordinate (-H/2 ≤ y ≤ H/2).")

# Execute the function to print the result.
display_final_equation()