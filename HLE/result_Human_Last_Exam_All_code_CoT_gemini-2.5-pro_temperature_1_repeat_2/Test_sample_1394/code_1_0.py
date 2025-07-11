import sympy

def solve_ode():
    """
    This function prints the general solution of the given differential equation.
    The solution was found by transforming the original equation into a Lagrange differential equation
    and then solving it.
    """
    x, y, k = sympy.symbols('x y k')
    
    # The general solution is y^2 = k*x^2 + k^2 + k + 9
    # where k is an arbitrary constant.
    
    lhs = y**2
    rhs = k*x**2 + k**2 + k + 9
    
    solution_eq = sympy.Eq(lhs, rhs)
    
    print("The general solution of the differential equation is:")
    print(solution_eq)
    
    # We can also print it out term by term as requested.
    print("\nIn the format: y**2 = k*x**2 + k**2 + k + 9")
    print(f"y**{2} = k*x**{2} + k**{2} + k + {9}")


solve_ode()