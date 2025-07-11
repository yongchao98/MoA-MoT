import sympy

def find_fixed_point_coupling():
    """
    Calculates and displays the leading order expression for the Wilson-Fisher fixed point
    coupling (u*) in phi^4 theory near D=4 dimensions.
    """
    # 1. Define the necessary symbolic variables.
    # u: the coupling constant
    # epsilon: the deviation from 4 dimensions (D = 4 - epsilon)
    # pi: the mathematical constant
    u, epsilon = sympy.symbols('u epsilon')
    pi = sympy.pi

    # 2. Define the one-loop beta function for phi^4 theory.
    # beta(u) = -epsilon*u + (3 * u**2) / (16 * pi**2)
    beta_func = -epsilon * u + (3 * u**2) / (16 * pi**2)
    
    print("The one-loop beta function is: β(u) =", sympy.pretty(beta_func, use_unicode=True))
    print("\n----------------------------------")

    # 3. A fixed point u* exists where β(u*) = 0. We solve this equation for u.
    print("Finding the fixed point u* by solving β(u*) = 0...")
    fixed_point_solutions = sympy.solve(beta_func, u)

    # The solutions include the trivial (u=0) and non-trivial (Wilson-Fisher) fixed points.
    # We are interested in the non-trivial one.
    u_star = None
    for sol in fixed_point_solutions:
        if sol != 0:
            u_star = sol
            break

    # 4. Display the final expression for the Wilson-Fisher fixed point u*.
    if u_star is not None:
        u_star_symbol = sympy.Symbol('u^*')
        final_equation = sympy.Eq(u_star_symbol, u_star)
        
        print("\nThe leading order expression for the Wilson-Fisher fixed point coupling is:")
        print(sympy.pretty(final_equation, use_unicode=True))

        print("\nBreaking down the components of the final equation:")
        # The expression is of the form (numerator * pi**2 / denominator) * epsilon
        numerator_coeff = sympy.fraction(u_star/epsilon)[0]/pi**2
        denominator_coeff = sympy.fraction(u_star/epsilon)[1]

        print(f"Numerical part of the coefficient: {numerator_coeff} / {denominator_coeff}")
        print(f"Term involving π: π^2")
        print(f"Dependence on ϵ: Proportional to ϵ")

    else:
        print("Could not find the non-trivial fixed point.")

if __name__ == '__main__':
    find_fixed_point_coupling()