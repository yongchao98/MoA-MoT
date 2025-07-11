def display_potential_distribution():
    """
    Prints the derived expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y) and its constituent parts.
    """

    # Define the symbols used in the equation as strings
    psi = "\u03C8"  # Greek letter psi
    y = "y"
    z1 = "z\u2081" # z with subscript 1
    beta = "\u03B2"   # Greek letter beta
    k = "k"
    H = "H"

    # Construct the components of the final equation
    amplitude = f"{z1}(1 + {beta}{k})"
    numerator = f"sinh(k(H/2 - {y}))"
    denominator = f"sinh(kH)"

    # Form the full expression string
    final_equation = f"{amplitude} * [{numerator} / {denominator}]"

    # Print the final result
    print("Based on the linearized Poisson-Boltzmann equation and the given boundary conditions,")
    print("the final expression for the Electrical Double-Layer (EDL) potential distribution is:")
    print(f"\n{psi}({y}) = {final_equation}\n")

    # Output each component of the final equation as requested
    print("The components of this equation are:")
    print(f"  - Amplitude Factor: {amplitude}")
    print(f"  - Numerator (describes spatial variation): {numerator}")
    print(f"  - Denominator (normalization factor): {denominator}")

# Execute the function to display the result
display_potential_distribution()