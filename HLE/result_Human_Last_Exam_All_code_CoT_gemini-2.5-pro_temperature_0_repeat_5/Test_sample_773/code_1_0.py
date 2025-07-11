import sympy

def solve_mass_calculation():
    """
    This function symbolically calculates the total mass based on the problem description.
    """
    # Define the symbols from the problem statement.
    # q_v is the order of the residual field of the local field K_hat.
    # q is the order of the constant field of the global function field.
    q_v = sympy.Symbol('q_v')
    q = sympy.Symbol('q')

    # The factor provided in the question.
    factor = (q_v * (q - 1)) / (q_v - 1)

    # Based on the theory of Tamagawa numbers and Haar measures, the volume
    # of the space of lattices, mu(X), is calculated to be 1/(q-1).
    mu_X = 1 / (q - 1)

    # The total mass is the product of the factor and the measure mu(X).
    total_mass = factor * mu_X

    # Simplify the expression to get the final result.
    simplified_mass = sympy.simplify(total_mass)

    # Print the final equation, showing each component.
    # We use sympy.pretty for a clear mathematical representation.
    factor_pretty = sympy.pretty(factor, use_unicode=False)
    mu_X_pretty = sympy.pretty(mu_X, use_unicode=False)
    result_pretty = sympy.pretty(simplified_mass, use_unicode=False)

    print("The final equation is the product of the given factor and the calculated measure mu(X):")
    print(f"({factor_pretty}) * ({mu_X_pretty}) = {result_pretty}")

if __name__ == '__main__':
    solve_mass_calculation()