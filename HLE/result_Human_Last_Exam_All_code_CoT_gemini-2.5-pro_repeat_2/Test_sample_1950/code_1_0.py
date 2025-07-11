import sympy

def solve_purification_protocol():
    """
    Calculates the product of the successful output fidelity and the success probability
    for the described GHZ state purification protocol.

    The calculation is performed symbolically using the sympy library.
    """
    # Define symbolic variables for the fidelities F1 and F2
    F1, F2 = sympy.symbols('F1 F2')

    # Define the coefficients for the input density matrices based on the problem description.
    # rho_GHZ(F1) = a1 * |GHZ><GHZ| + b1 * I_8
    # rho_Bell(F2) = a2 * |Phi+><Phi+| + b2 * I_4
    a1 = (8 * F1 - 1) / 7
    b1 = (1 - F1) / 7
    a2 = (4 * F2 - 1) / 3
    b2 = (1 - F2) / 3

    # The total product of output fidelity and success probability can be found by
    # decomposing the input state into four parts and summing their individual contributions.
    # The four parts arise from expanding the tensor product of the two mixed states.
    # The contribution from each part has been derived as:
    # 1. (pure GHZ, pure Bell): a1 * a2
    # 2. (pure GHZ, mixed Bell): a1 * b2
    # 3. (mixed GHZ, pure Bell): b1 * a2
    # 4. (mixed GHZ, mixed Bell): 2 * b1 * b2

    C1 = a1 * a2
    C2 = a1 * b2
    C3 = b1 * a2
    C4 = 2 * b1 * b2

    # The total product is the sum of these four contributions.
    total_product = C1 + C2 + C3 + C4

    # Simplify the final expression.
    simplified_product = sympy.simplify(total_product)

    print("The symbolic expression for the product of the output fidelity and success probability is:")
    sympy.pprint(simplified_product)

    # To fulfill the requirement of outputting each number in the final equation,
    # we can format the simplified expression into a clear string.
    # The simplified expression is (22*F1*F2 - F1 - F2 + 1) / 21
    print("\nThe final equation is:")
    print("(22 * F1 * F2 - 1 * F1 - 1 * F2 + 1) / 21")

    # Final answer in the specified format
    final_answer = str(simplified_product)
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_purification_protocol()