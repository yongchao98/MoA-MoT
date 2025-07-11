def solve_edl_potential():
    """
    This function constructs and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution, psi(y).
    """

    # Define the symbols used in the equation as strings for display.
    # Using Unicode for clarity where appropriate.
    psi_y = "ψ(y)"
    zeta_1 = "ζ₁"
    beta = "β"
    k = "k"
    H = "H"
    y = "y"
    one = "1"

    # Build the components of the final equation to show each part.
    # Term 1: The slip-dependent zeta potential at the bottom wall.
    term1 = f"{zeta_1}({one} + {beta}{k})"

    # Numerator of the hyperbolic ratio.
    numerator_arg = f"{k}({H} - {y})"
    numerator = f"sinh({numerator_arg})"

    # Denominator of the hyperbolic ratio.
    denominator_arg = f"{k}{H}"
    denominator = f"sinh({denominator_arg})"

    # Construct the final expression string.
    final_equation = f"{psi_y} = {term1} * ( {numerator} / {denominator} )"

    print("The final expression for the Electrical double-layer potential distribution is:")
    print(final_equation)

solve_edl_potential()

# The final answer in a simplified format.
# ψ(y) = ζ₁(1 + βk) * [sinh(k(H-y)) / sinh(kH)]
final_answer_string = "ψ(y) = ζ₁(1 + βk) * sinh(k(H-y)) / sinh(kH)"
print(f"\n<<< {final_answer_string} >>>")