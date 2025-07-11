import sympy as sp

def find_separatrix():
    """
    Finds the separatrix for the given system of differential equations by searching
    for an invariant manifold of the form d = a*u**2 + b*u + c.
    """
    # Define symbols
    t = sp.Symbol('t')
    u = sp.Function('u')(t)
    d = sp.Function('d')(t)
    a, b, c = sp.symbols('a b c')

    # The given system of differential equations
    u_prime = u**2 * (u - 1)
    # The original equation for d'(t)
    d_prime_orig = 2*d**2 + (-3*u + 5*u**2)*d - u**3 * (1 - u)

    # Postulate the form of the separatrix: d(t) = a*u(t)^2 + b*u(t) + c
    d_guess = a*u**2 + b*u + c

    # Calculate the derivative of the guess with respect to t using the chain rule:
    # d/dt [d_guess] = (d/du [d_guess]) * (du/dt)
    d_guess_prime = sp.diff(d_guess, u) * u_prime
    
    # Substitute the postulated form of d into the original equation for d'(t)
    d_prime_from_eq = d_prime_orig.subs(d, d_guess)

    # For d_guess to be an invariant manifold, the two expressions for d' must be equal
    # for all values of u. We create an equation where the difference is zero.
    equation = sp.expand(d_prime_from_eq - d_guess_prime)

    # The resulting expression is a polynomial in u. For this polynomial to be zero
    # for all u, each of its coefficients must be zero.
    poly = sp.Poly(equation, u)
    coeffs = poly.coeffs()

    # Solve the system of equations where each coefficient is set to zero.
    solutions = sp.solve(coeffs, (a, b, c))

    print("Found the following polynomial invariant manifolds of the form d = a*u^2 + b*u + c:")
    
    # Store the solution equations
    solution_equations = []
    u_sym = sp.Symbol('u')
    for sol in solutions:
        a_sol, b_sol, c_sol = sol[a], sol[b], sol[c]
        solution_equations.append(sp.Eq(sp.Symbol('d'), a_sol*u_sym**2 + b_sol*u_sym + c_sol))
        # Print the equation with all coefficients as requested
        print(f"d = ({a_sol})*u**2 + ({b_sol})*u + ({c_sol})")

    # The fixed points of the system are (0,0), (1,0) (unstable node), and (1,-1) (saddle).
    # Separatrices are the stable/unstable manifolds of saddle points.
    # We check which of our solutions passes through the saddle point (1, -1).
    saddle_point = {u_sym: 1, sp.Symbol('d'): -1}
    separatrix_eq = None
    final_coeffs = None
    
    for eq, sol in zip(solution_equations, solutions):
        # Check if eq.rhs - eq.lhs is zero at the saddle point
        if sp.simplify(eq.lhs - eq.rhs).subs(saddle_point) == 0:
            separatrix_eq = eq
            final_coeffs = sol

    if separatrix_eq:
        print("\nThe separatrix is the manifold that passes through the saddle point (1, -1).")
        a_sol, b_sol, c_sol = final_coeffs[a], final_coeffs[b], final_coeffs[c]
        print(f"The equation for the separatrix is:")
        print(f"d = ({a_sol})*u**2 + ({b_sol})*u + ({c_sol})")
    else:
        print("\nCould not identify a separatrix among the polynomial solutions.")

if __name__ == '__main__':
    find_separatrix()