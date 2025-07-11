import sympy

def solve_for_lower_bound():
    """
    Calculates the lower bound 'a' based on the analysis of the PDE.
    The lower bound is determined by the negative of the maximum value of the function f(u) = u(1-u)
    on the interval [0, 1].
    """
    # Define the variable u
    u = sympy.Symbol('u')

    # Define the function f(u) from the PDE's nonlinear term
    f_u = u * (1 - u)

    # To find the maximum value of f(u), we first find its derivative
    f_prime = sympy.diff(f_u, u)

    # Solve for the critical points by setting the derivative to zero
    critical_points = sympy.solve(f_prime, u)
    
    # We are interested in the interval [0, 1]
    u_max_val = critical_points[0]

    # Calculate the maximum value of f(u) at this point
    max_f_u = f_u.subs(u, u_max_val)

    # The lower bound 'a' is the negative of this maximum value
    a = -max_f_u

    print(f"The nonlinear term in the PDE is u*(1-u).")
    print(f"Let f(u) = u*(1-u) = u - u**2.")
    print(f"To find the maximum of f(u) for u in [0, 1], we compute the derivative: f'(u) = 1 - 2*u.")
    print(f"Setting f'(u) = 0 gives 1 - 2*u = 0, which means u = 1/2.")
    print(f"The maximum value of f(u) is at u = 1/2, which is f(1/2) = (1/2)*(1 - 1/2) = 1/4.")
    
    # Based on the theory of this type of equation, the lower bound is the negative of this value.
    numerator = -1
    denominator = 4
    
    print(f"The lower bound 'a' is the negative of this maximum value.")
    print(f"So, a = -1/4.")
    print(f"The final equation is:")
    print(f"a = {numerator} / {denominator}")
    print(f"a = {float(a)}")

solve_for_lower_bound()