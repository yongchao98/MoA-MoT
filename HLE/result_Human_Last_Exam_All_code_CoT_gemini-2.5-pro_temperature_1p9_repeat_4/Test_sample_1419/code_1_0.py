import sympy as sp

def solve_fixed_point():
    """
    Calculates the leading order expression for the Wilson-Fisher fixed point
    coupling u* in phi^4 theory near d=4 dimensions.
    """
    # Define the symbols for the coupling 'u' and the small parameter 'epsilon'
    u, epsilon, pi = sp.symbols('u epsilon pi')

    # The one-loop beta function for phi^4 theory in d = 4 - epsilon dimensions.
    # The convention used is for a Lagrangian with an interaction term (u/4!) * phi^4.
    # The coefficient is 3 / (16 * pi^2) for a single-component scalar field (N=1).
    coefficient = 3 / (16 * pi**2)
    beta_u = -epsilon * u + coefficient * u**2

    print("The one-loop beta function for the coupling u in ϕ⁴ theory is:")
    # Using symbols for a cleaner printout
    print(f"β(u) = -ε*u + (3 / (16*π²)) * u²\n")

    # A fixed point u* is found where the beta function is zero: β(u*) = 0
    # We are interested in the non-trivial fixed point (u* != 0).
    # The equation is: -ε*u* + (3 / (16*π²)) * (u*)² = 0
    # For u* != 0, we can divide by u*: -ε + (3 / (16*π²)) * u* = 0
    
    fixed_point_eq = sp.Eq(-epsilon + coefficient * u, 0)

    print("The equation for the non-trivial fixed point u* is:")
    # To avoid floating point numbers in the printout, build the expression with Integers
    u_sym = sp.Symbol('u*')
    eq_to_print = sp.Eq(-epsilon + (sp.Integer(3) / (sp.Integer(16) * pi**2)) * u_sym, 0)
    print(f"{sp.pretty(eq_to_print, use_unicode=True)}\n")

    # Solve the equation for u to find the expression for the fixed point u*
    u_star_solution = sp.solve(fixed_point_eq, u)
    
    # The result from solve is a list, so we take the first element
    u_star = u_star_solution[0]

    # Now, let's explicitly print the final expression for u* with all its numerical parts.
    num_coeff = 16
    den_coeff = 3
    print("Solving for u*, we get the leading order expression for the Wilson-Fisher fixed point:")
    print(f"u* = ({num_coeff} * π² * ε) / {den_coeff}")

if __name__ == '__main__':
    solve_fixed_point()
