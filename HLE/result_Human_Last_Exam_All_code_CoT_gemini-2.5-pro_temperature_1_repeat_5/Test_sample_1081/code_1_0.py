def print_e8_tori_formula():
    """
    Constructs and prints the formula for the number of F_q-rational maximal tori
    in a reductive group of type E_8.
    """
    # Parameters for the E_8 group structure
    rank = 8
    num_positive_roots = 120
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]

    # Construct the numerator of the formula
    numerator_parts = [f"(q^{d} - 1)" for d in degrees]
    numerator_str = f"q^{num_positive_roots} * {' * '.join(numerator_parts)}"

    # Construct the denominator of the formula
    denominator_str = f"(q - 1)^{rank}"

    # Final formula
    final_formula = f"({numerator_str}) / ({denominator_str})"

    print("The exact number of F_q-rational maximal tori of a group of type E_8 is given by the formula:")
    print(final_formula)
    print("\nAs requested, here are the specific numbers that appear in the equation, based on the E_8 root system:")
    print(f"- The exponent in the denominator is the rank of E_8: {rank}")
    print(f"- The exponent of the leading q term in the numerator is the number of positive roots of E_8: {num_positive_roots}")
    print(f"- The exponents in the factors (q^d - 1) are the degrees of the fundamental invariants of E_8: {degrees}")

print_e8_tori_formula()