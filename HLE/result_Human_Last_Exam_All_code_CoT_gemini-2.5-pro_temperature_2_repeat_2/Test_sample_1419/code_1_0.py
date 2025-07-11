import sympy

def find_wilson_fisher_fixed_point():
    """
    This function symbolically derives the leading order expression for the
    Wilson-Fisher fixed point coupling in phi^4 theory near d=4 dimensions.
    """
    # Define the symbolic variables. We use Greek letters for readability.
    # u_star will represent the fixed-point coupling constant
    # epsilon represents the deviation from 4 dimensions (d = 4 - epsilon)
    # pi is the mathematical constant
    u_star, epsilon, pi = sympy.symbols('u* ε π')

    # The one-loop beta function for the dimensionless coupling 'u' in phi^4 theory is:
    # β(u) = -ε*u + (3*u^2)/(16*π^2)
    # This equation describes how the coupling 'u' changes with the energy scale.
    beta_function = -epsilon * u_star + (3 * u_star**2) / (16 * pi**2)

    # A fixed point is a point where the coupling does not change with scale, i.e., β(u*) = 0.
    # We solve this equation to find the fixed point coupling u*.
    # We are interested in the non-trivial solution (where u* is not zero),
    # which is the Wilson-Fisher fixed point.
    fixed_point_solutions = sympy.solve(beta_function, u_star)

    # sympy.solve returns a list of solutions. For β(u*)=0, the solutions are
    # [0, (16*π²*ε)/3]. The first solution is the trivial Gaussian fixed point.
    # The second is the non-trivial Wilson-Fisher fixed point.
    wilson_fisher_fixed_point = fixed_point_solutions[1]

    # Now we have the symbolic expression for the fixed point.
    # The prompt requires us to output each number in the final equation.
    # We can extract the numerical parts of the coefficient from the expression.
    # The expression is of the form (numerator/denominator) * π² * ε
    coeff = wilson_fisher_fixed_point / (pi**2 * epsilon)
    numerator, denominator = coeff.as_numer_denom()

    print("The leading order expression for the fixed point coupling u* is derived by setting the one-loop β-function to zero.")
    print("The β-function is: β(u) = -ε*u + 3*u²/(16*π²)")
    print("Setting β(u*) = 0 and solving for the non-trivial fixed point gives the following equation:")
    
    # Print the final equation with numbers separated as requested
    print("\nu* = (", int(numerator), " * π² /", int(denominator), ") * ε")

if __name__ == "__main__":
    find_wilson_fisher_fixed_point()