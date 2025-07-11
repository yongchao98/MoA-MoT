import sympy

def solve_ladder_cutoff_frequency():
    """
    Calculates the cutoff frequency at node a0 in the given infinite ladder network
    by finding the Thevenin equivalent resistance symbolically using sympy.
    """
    # Define symbolic variables for resistance, capacitance, and other helper variables
    r, C = sympy.symbols('r C', positive=True, real=True)
    A, B = sympy.symbols('A B') # Coefficients for the general solution

    # The decaying mode eigenvalue for the recurrence relations is lambda = 2 - sqrt(3)
    lambda_val = 2 - sympy.sqrt(3)

    # Based on the analysis of the recurrence relations and boundary conditions,
    # a system of equations for the coefficients A and B is derived.
    # Equation 1 from KCL at node c1: A*(4 - lambda) + B = 0
    # Equation 2 from KCL at node d1: A*(lambda - 4) + B = r
    
    eq1 = sympy.Eq(A * (4 - lambda_val) + B, 0)
    eq2 = sympy.Eq(A * (lambda_val - 4) + B, r)
    
    # Solve this system for A and B
    solution = sympy.solve([eq1, eq2], (A, B))
    A_sol = solution[A]
    B_sol = solution[B]

    # The Thevenin resistance R_th is the voltage at node d0 (=a0) for a 1A test current.
    # v_d0 = v_d1 + r, where v_d1 is the voltage at node d1.
    # For n=1, the voltage v_d1 is given by -A + B from the general solution form.
    v_d1 = -A_sol + B_sol
    R_th = v_d1 + r
    
    # Simplify the resulting expression for R_th
    R_th_simplified = sympy.simplify(R_th)

    # The cutoff angular frequency is omega_c = 1 / (R_th * C)
    omega_c = 1 / (R_th_simplified * C)
    omega_c_simplified = sympy.simplify(omega_c)

    # Deconstruct the final expression to print it clearly
    num, den = sympy.fraction(omega_c_simplified)
    # The denominator is proportional to r*C. Let's extract the constant factor.
    den_factor = sympy.simplify(den / (r * C))
    
    # The factor is (5 - sqrt(3)). We can extract the integer parts.
    c1 = 5
    c2 = 3
    
    print(f"The Thevenin resistance R_th seen from node a0 is: r * ({c1} - sqrt({c2})) / 2")
    print("\nThe cutoff angular frequency omega_c = 1 / (R_th * C) is therefore:")
    print(f"omega_c = {num} / (({c1} - sqrt({c2})) * r * C)")


solve_ladder_cutoff_frequency()