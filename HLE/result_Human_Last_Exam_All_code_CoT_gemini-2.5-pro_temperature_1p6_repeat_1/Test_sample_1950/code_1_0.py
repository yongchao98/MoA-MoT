import sympy

def solve_ghz_purification():
    """
    Calculates the product of the success probability and the output fidelity
    for the described GHZ state purification protocol.
    The result is a symbolic expression in terms of the initial fidelities F1 and F2.
    """
    # Define F1 and F2 as symbolic variables
    F1, F2 = sympy.symbols('F1 F2')

    # Define the probabilities p1 and p2 based on the input state definitions
    # p1 is the weight of the pure GHZ state in the input mixture
    p1 = (8 * F1 - 1) / 7
    # p2 is the weight of the pure Bell state in the input mixture
    p2 = (4 * F2 - 1) / 3

    # The product of success probability and output fidelity can be calculated by
    # summing the contributions from the four components of the tensor product of the input density matrices.
    # Let G = |GHZ><GHZ| and B = |Bell><Bell|.
    # The input state is (p1*G + (1-p1)*I/8) @ (p2*B + (1-p2)*I/4)
    #
    # 1. Contribution from the (G @ B) term with coefficient p1*p2.
    #    The ideal state component remains ideal. The fidelity contribution is 1.
    term1 = p1 * p2 * 1

    # 2. Contribution from the (G @ I/4) term with coefficient p1*(1-p2).
    #    The output state is partially decohered, yielding a fidelity contribution of 1/4.
    term2 = p1 * (1 - p2) * sympy.Rational(1, 4)

    # 3. Contribution from the (I/8 @ B) term with coefficient (1-p1)*p2.
    #    The operations depolarize the Bell state, resulting in a fidelity contribution of 1/16.
    term3 = (1 - p1) * p2 * sympy.Rational(1, 16)

    # 4. Contribution from the (I/8 @ I/4) term with coefficient (1-p1)*(1-p2).
    #    The maximally mixed state results in a fidelity contribution of 1/16.
    term4 = (1 - p1) * (1 - p2) * sympy.Rational(1, 16)

    # The total product is the sum of these four terms.
    product_fidelity_prob = term1 + term2 + term3 + term4

    # Simplify the final expression
    simplified_product = sympy.simplify(product_fidelity_prob)

    # To satisfy the "output each number in the final equation" requirement,
    # we manually construct the output string from the known simplified expression.
    # The simplified expression is (16*F1*F2 - F1 - 2*F2 + 1) / 14
    
    print("Let F1 be the initial fidelity of the 3-qubit GHZ state.")
    print("Let F2 be the initial fidelity of the 2-qubit Bell state.")
    print("The product of the successful output fidelity and the success probability is given by the expression:")
    print("(16 * F1 * F2 - 1 * F1 - 2 * F2 + 1) / 14")


solve_ghz_purification()