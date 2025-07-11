import sympy

def solve_purification_product():
    """
    This function calculates and prints the symbolic expression for the product of
    the success probability and the output fidelity of the GHZ purification protocol.
    """

    # Step 1: Define symbolic variables for the input fidelities F1 and F2.
    F1, F2 = sympy.symbols('F1 F2')

    # Step 2: Define the mixture parameters p1 and p2 based on the definitions
    # of the input Werner states.
    # p1 corresponds to the 3-qubit GHZ state with fidelity F1.
    p1 = (8 * F1 - 1) / 7

    # p2 corresponds to the 2-qubit Bell state with fidelity F2.
    p2 = (4 * F2 - 1) / 3

    # Step 3: Use the derived formula for the product of the success probability
    # and the output fidelity in terms of p1 and p2.
    product_in_p = (1 + 3 * p1 + p2 + 11 * p1 * p2) / 16

    # Step 4: Substitute the expressions for p1 and p2 into the formula and
    # simplify it to obtain the final expression in terms of F1 and F2.
    final_formula = sympy.simplify(product_in_p)

    # Step 5: Print the final result. The simplified expression clearly shows
    # all the numbers involved in the final equation as requested.
    print("The product of the successful output fidelity and the success probability is:")
    print(final_formula)

    # To explicitly output each number as requested by the prompt,
    # we can extract the coefficients of the polynomial.
    num, den = final_formula.as_numer_denom()
    num_poly = sympy.Poly(num, F1, F2)
    c_const = num_poly.coeff_monomial(1)
    c_f1 = num_poly.coeff_monomial(F1)
    c_f2 = num_poly.coeff_monomial(F2)
    c_f1f2 = num_poly.coeff_monomial(F1*F2)

    print("\nThe equation has the form (c0 + c1*F1 + c2*F2 + c12*F1*F2) / d")
    print(f"The numbers in the equation are:")
    print(f"c0 = {c_const}")
    print(f"c1 = {c_f1}")
    print(f"c2 = {c_f2}")
    print(f"c12 = {c_f1f2}")
    print(f"d = {den}")


solve_purification_product()