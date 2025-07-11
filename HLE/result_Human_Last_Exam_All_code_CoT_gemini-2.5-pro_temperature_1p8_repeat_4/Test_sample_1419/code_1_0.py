import sympy

def find_fixed_point_coupling():
    """
    Calculates and prints the leading order expression for the Wilson-Fisher
    fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    # 1. Define the symbolic variables used in the equation.
    # u_star represents the fixed point coupling u*.
    # epsilon represents the deviation from 4 dimensions (d = 4 - epsilon).
    # pi is the mathematical constant.
    u_star = sympy.Symbol('u^*')
    epsilon = sympy.Symbol('epsilon')
    pi = sympy.pi

    # 2. Define the one-loop beta function for phi^4 theory.
    # beta(u) = -epsilon*u + (3*u^2)/(16*pi^2)
    # At the fixed point u*, beta(u*) = 0.
    beta_function_at_fixed_point = -epsilon * u_star + (3 * u_star**2) / (16 * pi**2)
    
    # 3. Create the equation to solve: beta(u*) = 0
    equation = sympy.Eq(beta_function_at_fixed_point, 0)
    
    # 4. Solve the equation for u_star.
    # The equation has two solutions:
    # a) u* = 0 (the trivial Gaussian fixed point)
    # b) The non-trivial Wilson-Fisher fixed point.
    solutions = sympy.solve(equation, u_star)
    
    # 5. Extract the non-trivial solution. sympy.solve returns a list.
    # The first solution is 0, the second is the one we want.
    wilson_fisher_fixed_point = solutions[1]
    
    # 6. Print the final expression for the fixed point coupling.
    # The result is u* = (16 * pi**2 / 3) * epsilon.
    # We will print the components of this equation as requested.
    
    coefficient = wilson_fisher_fixed_point / epsilon
    num, den = coefficient.as_numer_denom()

    # The problem asks to output each number in the final equation.
    print("The equation for the fixed point coupling u^* is:")
    print(f"{u_star} = ({num.evalf()} / {den.evalf()}) * {pi**2} * {epsilon}")
    
    # Let's print it more symbolically and clearly.
    # The components are 'u^*', '=', '(', '16', '*', 'pi**2', '/', '3', ')', '*', 'epsilon'
    print("\nSymbolic expression:")
    print(str(u_star), "=", "(", str(sympy.numer(coefficient)), "*", "pi**2", "/", str(sympy.denom(coefficient)), ")", "*", str(epsilon))

if __name__ == '__main__':
    find_fixed_point_coupling()