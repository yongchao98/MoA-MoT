import sympy

def solve_dessin_equations():
    """
    Solves the system of equations derived from the properties of the simple dessin.
    The system is:
    1) d = 8
    2) d = (N_r - 1) * 2 * m
    3) d = N_r * m
    where d, N_r, m are positive integers and m is even.
    We are looking for the value of N_r.
    """
    d = 8
    N_r = sympy.Symbol('N_r', integer=True, positive=True)
    m = sympy.Symbol('m', integer=True, positive=True)
    
    # From d = N_r * m
    eq1 = sympy.Eq(d, N_r * m)
    
    # From d = (N_r - 1) * 2 * m
    eq2 = sympy.Eq(d, (N_r - 1) * 2 * m)
    
    # We can solve this system.
    # From eq1: m = d / N_r = 8 / N_r
    # Substitute into eq2: 8 = (N_r - 1) * 2 * (8 / N_r)
    # 8 * N_r = (N_r - 1) * 16
    # 8 * N_r = 16 * N_r - 16
    # 8 * N_r = 16
    # N_r = 2
    
    # Now find m using N_r = 2
    # 8 = 2 * m  => m = 4
    
    # Let's use sympy to solve it more formally
    solutions = sympy.solve([eq1, eq2], (N_r, m))
    
    # The solver might return multiple solutions, let's filter for our constraints.
    final_N_r = None
    final_m = None
    
    if solutions:
        sol_N_r = solutions[N_r]
        sol_m = solutions[m]
        # Check constraint: m must be a positive even integer.
        if sol_m > 0 and sol_m % 2 == 0:
            final_N_r = sol_N_r
            final_m = sol_m

    # Output the results in a readable format, including the equation.
    if final_N_r is not None:
        print(f"Solving the system of equations for the degree d={d}:")
        print(f"Equation from q-vertices: {d} = ({final_N_r} - 1) * 2 * {final_m}")
        print(f"Equation from r-vertices: {d} = {final_N_r} * {final_m}")
        
        # Check if equations hold
        q_check = (final_N_r - 1) * 2 * final_m
        r_check = final_N_r * final_m
        print(f"Verification: {q_check} = {d} and {r_check} = {d}")

        print("\nThe determined values are:")
        print(f"Number of r-vertices in ]0, 1[ (N_r) = {final_N_r}")
        print(f"Valency parameter (m) = {final_m}")

        print("\nThe maximum number of vertices labelled r within ]0, 1[ is:")
        print(final_N_r)
    else:
        print("No solution found that satisfies all constraints.")

solve_dessin_equations()