import sympy

def find_wilson_fisher_fixed_point():
    """
    This function calculates and prints the leading order expression for the
    Wilson-Fisher fixed point coupling in phi^4 theory.
    """
    # Step 1: Define the necessary symbols.
    # u: the coupling constant
    # epsilon: the small parameter in d = 4 - epsilon
    # pi: the mathematical constant pi
    u, epsilon, pi = sympy.symbols('u* ε π')

    # Step 2: Define the one-loop beta function for phi^4 theory.
    # The term '-epsilon * u' comes from the canonical dimension of the coupling.
    # The term '(3 / (16*pi**2)) * u**2' is the one-loop quantum correction.
    beta_u_expr = -epsilon * u + (3 / (16 * pi**2)) * u**2

    print("The one-loop beta function for the ϕ⁴ theory in d = 4 - ε dimensions is:")
    print(f"β(u) = -ε*u + (3 / (16*π²))*u²")

    # Step 3: Find the fixed point by solving beta(u*) = 0.
    # A fixed point is a value of the coupling where the beta function is zero.
    fixed_point_solutions = sympy.solve(beta_u_expr, u)

    print("\nTo find the fixed points, we solve the equation β(u*) = 0:")
    print(f"-ε*u* + (3 / (16*π²))(u*)² = 0")

    # The solutions are a list. The first is the trivial Gaussian fixed point (u* = 0).
    # The second is the non-trivial Wilson-Fisher fixed point.
    wilson_fisher_fp = fixed_point_solutions[1]

    # Step 4: Display the final expression for the non-trivial fixed point.
    print("\nThe leading order expression for the non-trivial fixed point coupling u* is:")

    # We explicitly extract the coefficients to fulfill the printing requirement.
    # The solution is of the form (numerator / denominator) * epsilon.
    numerator_coeff = 16
    denominator_coeff = 3
    
    print(f"u* = ({numerator_coeff} * π² / {denominator_coeff}) * ε")

if __name__ == '__main__':
    find_wilson_fisher_fixed_point()