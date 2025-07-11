import sympy

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # Define symbols and the function
    n = sympy.symbols('n', integer=True)
    y = sympy.Function('y')

    # The equation is 8y[n] - 6y[n-1] + y[n-2] = 1
    # For rsolve, we write it as an expression that equals 0
    # 8y(n) - 6y(n-1) + y(n-2) - 1 = 0
    equation = 8*y(n) - 6*y(n-1) + y(n-2) - 1

    # The initial conditions are y[0] = 1 and y[-1] = 2.
    # rsolve doesn't handle negative indices well. We can find y[1] to use as an initial condition instead.
    # At n=1: 8*y[1] - 6*y[0] + y[-1] = 1
    #          8*y[1] - 6*(1) + 2 = 1
    #          8*y[1] - 4 = 1
    #          8*y[1] = 5  => y[1] = 5/8
    initial_conditions = {y(0): 1, y(1): sympy.S(5)/8}

    # Solve the recurrence equation
    solution = sympy.rsolve(equation, y(n), initial_conditions)

    # The solution is of the form y[n] = A*(B)**n + C*(D)**n + E
    # We need to parse the symbolic solution to find the constants.

    # The constant term E is the part of the solution that does not depend on n.
    E = solution.subs(n, sympy.oo)

    # The homogeneous part is the rest of the solution.
    hom_solution = solution - E

    # Helper function to parse terms of the form Coeff * Base**n
    def parse_term(term, n_sym):
        coeff = term.subs(n_sym, 0)
        base_n_part = term / coeff
        base = base_n_part.subs(n_sym, 1)
        return coeff, base

    # Get the two terms from the homogeneous solution
    hom_terms = hom_solution.as_ordered_terms()
    
    term1_coeff, term1_base = parse_term(hom_terms[0], n)
    term2_coeff, term2_base = parse_term(hom_terms[1], n)

    # Assign A, B, C, D based on the convention B > D
    if term1_base > term2_base:
        A, B = term1_coeff, term1_base
        C, D = term2_coeff, term2_base
    else:
        A, B = term2_coeff, term2_base
        C, D = term1_coeff, term1_base

    # Print the equation with the found constants
    print("The closed-form solution is:")
    print(f"y[n] = {A} * ({B})^n + {C} * ({D})^n + {E}\n")

    # Print the individual constants
    print(f"The constants are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}\n")
    
    # Calculate the final expression
    result = E/A + (D*C)/B

    # Print the result
    print(f"The value of E/A + (D*C)/B is:")
    print(f"({E})/({A}) + (({D})*({C}))/({B}) = {result.evalf()}")

solve_difference_equation()
<<<0.75>>>