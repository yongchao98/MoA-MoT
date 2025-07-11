import sympy

def solve_fixed_point_coupling():
    """
    This script calculates and displays the leading order expression for the
    Wilson-Fisher fixed point coupling in phi^4 theory.
    """
    # Define symbols for the equation.
    u_star, epsilon, pi = sympy.symbols('u^* epsilon pi')

    # The fixed point equation is beta(u^*) = 0.
    # At one-loop order, the beta function for the coupling u is:
    # beta(u) = -epsilon * u + (3 * u**2) / (16 * pi**2)
    # Setting this to zero for u = u_star gives the equation to solve.
    fixed_point_eq = sympy.Eq(-epsilon * u_star + (3 * u_star**2) / (16 * pi**2), 0)

    # Solve the equation for u_star. The solver will return both the trivial
    # and non-trivial solutions.
    solutions = sympy.solve(fixed_point_eq, u_star)

    # We are interested in the non-trivial (Wilson-Fisher) fixed point.
    # solutions list is [0, 16*pi**2*epsilon/3]
    nontrivial_solution = solutions[1]
    
    # Extract the numerical components of the expression for printing.
    # The expression is u^* = (16 * pi**2 / 3) * epsilon.
    numerator_val = 16
    symbol_str = "pi"
    exponent_val = 2
    denominator_val = 3
    
    print("In phi^4 theory near four dimensions, the non-trivial fixed point u^* is found by solving the equation beta(u^*) = 0.")
    print("\nThe one-loop beta function leads to the equation:")
    print(f"-epsilon * u^* + (3 * u^*^2) / (16 * pi^2) = 0")
    
    print("\nSolving for u^* yields the leading order expression for the Wilson-Fisher fixed point.")
    print("The final equation is composed of the following numbers and symbols:")
    print(f"  - A numerator coefficient: {numerator_val}")
    print(f"  - A symbolic constant: {symbol_str}")
    print(f"  - An exponent for the constant: {exponent_val}")
    print(f"  - A denominator: {denominator_val}")

    print("\nThe final expression is:")
    print(f"u^* = ({numerator_val} * {symbol_str}**{exponent_val} / {denominator_val}) * epsilon")

if __name__ == '__main__':
    solve_fixed_point_coupling()