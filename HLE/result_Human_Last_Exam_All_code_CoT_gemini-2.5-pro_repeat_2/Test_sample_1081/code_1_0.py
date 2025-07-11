def solve_tori_number():
    """
    This function calculates and prints the formula for the number of
    F_q-rational maximal tori of a reductive group of type E_8.
    """
    # The degrees of the fundamental invariants for the Weyl group of type E_8.
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]

    # We construct the formula for the number of F_q-rational maximal tori.
    # The formula is the square of the product of (q^{d_i} - 1) for each degree d_i.

    # Let's build the product term string
    product_terms = []
    for d in degrees:
        product_terms.append(f"(q^{d} - 1)")

    product_string = " * ".join(product_terms)

    # The final equation for the number of tori, let's call it N
    final_equation = f"N = ({product_string})^2"

    print("The exact number (N) of F_q-rational maximal tori for a group of type E_8 over F_q is given by the equation:")
    print(final_equation)

    print("\nThe numbers appearing in the exponents of the equation are the degrees of the fundamental invariants of the Weyl group W(E_8):")
    # Now, print each number as requested
    for d in degrees:
        print(d)

solve_tori_number()