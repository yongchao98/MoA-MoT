def display_potential_distribution():
    """
    This function prints the derived mathematical expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y).
    """

    # Define the symbols as strings for readable output
    z_1 = "z_1"
    beta = "β"
    k = "k"
    H = "H"
    y = "y"

    # The problem specifies the slip-dependant zeta potential z_a1
    # which acts as a coefficient in the final expression.
    # z_a1 = z_1*(1 + beta*k)
    coefficient = f"{z_1}*(1 + {beta}*{k})"

    # The numerator of the hyperbolic fraction
    # sinh(k*(H/2 - y))
    numerator = f"sinh({k}*({H}/2 - {y}))"

    # The denominator of the hyperbolic fraction
    # sinh(k*H)
    denominator = f"sinh({k}*{H})"

    # Construct the final expression for psi(y)
    final_equation = f"ψ({y}) = {coefficient} * [{numerator} / {denominator}]"

    # Print the final result clearly
    print("The final expression for the Electrical double-layer potential distribution is:")
    print(final_equation)

# Execute the function to display the result
display_potential_distribution()