def display_potential_distribution_expression():
    """
    Prints the derived expression for the electrical double-layer (EDL) potential distribution.
    """
    expression_parts = {
        "psi(y)": "The electrical potential at a vertical distance 'y' from the bottom plate.",
        "zeta_1": "The base zeta potential at the bottom surface (j=1) before slip consideration.",
        "beta": "The slip length parameter.",
        "k": "The Debye-Hückel parameter.",
        "H": "The height of the microchannel.",
        "y": "The vertical coordinate, ranging from 0 (bottom plate) to H (top plate).",
        "sinh": "The hyperbolic sine function."
    }

    print("The final expression for the Electrical Double-Layer (EDL) potential distribution psi(y) is:")
    print("========================================================================================\n")
    print("psi(y) = zeta_a1 * [sinh(k * (H - y)) / sinh(k * H)]")
    print("\nwhere the slip-dependant zeta potential at the bottom surface, zeta_a1, is given by:")
    print("zeta_a1 = zeta_1 * (1 + beta * k)")
    print("\nSo the full expression is:\n")
    # This prints the final equation with all components, as requested.
    print("psi(y) = zeta_1 * (1 + beta * k) * [sinh(k * (H - y)) / sinh(k * H)]")
    print("\n========================================================================================")
    print("Description of the variables in the equation:")
    for var, desc in expression_parts.items():
        print(f"- {var}: {desc}")
    print("\nNote: This expression is derived under the condition that the zeta potential at the top surface, zeta_2, is zero.")


if __name__ == "__main__":
    display_potential_distribution_expression()