import sympy

def solve_purification_protocol():
    """
    Calculates the product of success probability and output fidelity
    for the described GHZ purification protocol.
    """
    # Define F1 and F2 as symbolic variables
    F1, F2 = sympy.symbols('F1 F2')

    # Define the coefficients for the input states as per the problem description.
    # For rho_GHZ(F1):
    alpha1 = (8 * F1 - 1) / 7
    beta1 = (1 - F1) / 7
    # For rho_Bell(F2):
    alpha2 = (4 * F2 - 1) / 3
    beta2 = (1 - F2) / 3

    # As derived from the step-by-step analysis, the product of the successful
    # output fidelity and the success probability is given by the sum of four terms
    # corresponding to the different parts of the input density matrix.
    # Contribution from (GHZ, Bell) part: alpha1 * alpha2
    # Contribution from (GHZ, I) part: alpha1 * beta2
    # Contribution from (I, Bell) part: beta1 * alpha2
    # Contribution from (I, I) part: 2 * beta1 * beta2
    product = alpha1 * alpha2 + alpha1 * beta2 + beta1 * alpha2 + 2 * beta1 * beta2

    # Simplify the expression for the final answer.
    # This can also be expressed more concisely as: F1*F2 + beta1*beta2
    simplified_product = sympy.simplify(product)

    # The problem asks to output each number in the final equation.
    # Printing the simplified symbolic expression achieves this by showing
    # all the numerical coefficients and terms involving F1 and F2.
    print(simplified_product)

if __name__ == '__main__':
    solve_purification_protocol()