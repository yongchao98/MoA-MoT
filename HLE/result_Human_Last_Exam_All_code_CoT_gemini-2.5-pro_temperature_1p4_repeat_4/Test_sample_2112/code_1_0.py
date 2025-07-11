import sympy as sp

def solve_for_r0():
    """
    Solves for the radial distance r0 by finding the roots of the derived algebraic equations
    and filtering for the one that satisfies r0 > 15.
    """
    r0 = sp.symbols('r0')
    
    # Define the functions f(r) and g(r)
    f_r0 = (4*r0 + 37) / (3 - r0)
    g_r0 = (3*r0 - 37) / (r0 + 4)
    
    # The conditions are |f(r0)| = 1/sqrt(2) or |g(r0)| = 1/sqrt(2)
    c = 1 / sp.sqrt(2)
    
    equations = [
        sp.Eq(f_r0, c),
        sp.Eq(f_r0, -c),
        sp.Eq(g_r0, c),
        sp.Eq(g_r0, -c)
    ]
    
    potential_solutions = []
    for eq in equations:
        potential_solutions.extend(sp.solve(eq, r0))
        
    # Find the unique solution greater than 15
    valid_solution = None
    for sol in potential_solutions:
        # Check if the solution is real and greater than 15
        if sol.is_real and sol.evalf() > 15:
            valid_solution = sol
            break
            
    if valid_solution:
        # Rationalize the expression to get a cleaner form
        # The valid solution comes from g(r0) = 1/sqrt(2)
        # r0 = (4*c+37)/(3-c) = (4/sqrt(2)+37)/(3-1/sqrt(2)) = (4+37*sqrt(2))/(3*sqrt(2)-1)
        # To rationalize, multiply numerator and denominator by (3*sqrt(2)+1)
        # Numerator: (4+37*sqrt(2))*(3*sqrt(2)+1) = 12*sqrt(2)+4+37*2*3+37*sqrt(2) = 222+4+49*sqrt(2) = 226+49*sqrt(2)
        # Denominator: (3*sqrt(2)-1)*(3*sqrt(2)+1) = 18-1 = 17
        num_1 = 226
        num_2 = 49
        num_3 = 2
        num_4 = 17
        final_eq = f"({num_1} + {num_2} * sqrt({num_3})) / {num_4}"

        print("The final expression for r0 is:")
        print(f"r0 = {final_eq}")
        print("\nThe numbers in the final equation are:")
        print(f"Number 1: {num_1}")
        print(f"Number 2: {num_2}")
        print(f"Number 3: {num_3} (under the square root)")
        print(f"Number 4: {num_4}")
        
        # Calculate the numerical value
        numerical_value = valid_solution.evalf()
        print(f"\nThe numerical value of r0 is: {numerical_value}")
    else:
        print("No solution found for r0 > 15.")

solve_for_r0()