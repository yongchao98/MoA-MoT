import sympy as sp

def solve_potential():
    """
    This function defines the symbolic variables and constructs the expression
    for the electric potential in the region 0 <= y <= a based on the derivation.
    It then prints the final formula.
    """
    # Define the symbolic variables used in the problem
    sigma_0, k, x, y, a, b, epsilon_1, epsilon_2 = sp.symbols(
        'sigma_0 k x y a b epsilon_1 epsilon_2', real=True, positive=True
    )

    # Based on the step-by-step derivation, we construct the final expression
    # for the potential in the region 0 < y < a.

    # Denominator term
    denominator = k * (epsilon_2 * sp.cosh(k * a) * sp.sinh(k * b) +
                       epsilon_1 * sp.sinh(k * a) * sp.cosh(k * b))

    # Numerator term
    numerator = -sigma_0 * sp.sinh(k * b) * sp.sinh(k * (y - a)) * sp.sin(k * x)

    # The full expression for the potential Phi in the specified region
    potential_phi = numerator / denominator

    print("The electric potential Phi(x, y) in the region 0 <= y <= a is:")
    # The sp.pretty() function provides a more readable output for the equation
    final_equation_str = f"Phi(x, y) = {sp.pretty(potential_phi, use_unicode=False)}"
    
    # Printing the final equation with all its symbolic components.
    print(final_equation_str)
    print("\nThis result matches option A.")


if __name__ == '__main__':
    solve_potential()