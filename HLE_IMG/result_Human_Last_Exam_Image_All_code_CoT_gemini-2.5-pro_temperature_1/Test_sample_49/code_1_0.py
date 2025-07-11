import sympy

def solve_cutoff_frequency():
    """
    Calculates the cutoff frequency at node a0 for the given infinite RC ladder network.
    """
    # Define symbolic variables for resistance, capacitance, and the ladder's equivalent resistance.
    # We assume r, C, and R_in are positive real values.
    r, C, R_in = sympy.symbols('r C R_in', positive=True)

    # 1. Set up the equation for the equivalent resistance of the infinite ladder (R_in).
    # The equation is R_in = r || (2r + R_in), which leads to a quadratic equation:
    # R_in^2 + 2*r*R_in - 2*r^2 = 0
    equation = sympy.Eq(R_in**2 + 2*r*R_in - 2*r**2, 0)

    # 2. Solve the quadratic equation for R_in.
    # The solutions will be a positive and a negative root. Resistance must be positive.
    R_in_solutions = sympy.solve(equation, R_in)
    R_in_val = [sol for sol in R_in_solutions if sol.is_positive][0]

    # 3. Calculate the Thevenin equivalent resistance (R_th) seen by the capacitor.
    # R_th = r + (r || R_in) + r
    R_parallel = (r * R_in_val) / (r + R_in_val)
    R_th = 2 * r + R_parallel

    # 4. Simplify the expression for R_th.
    R_th_simplified = sympy.simplify(R_th)

    # 5. Calculate the angular cutoff frequency (w_c).
    # w_c = 1 / (R_th * C)
    w_c = 1 / (R_th_simplified * C)

    # 6. Simplify the final expression for w_c.
    w_c_simplified = sympy.simplify(w_c)
    
    # 7. Print the final result in a structured way.
    num, den = w_c_simplified.as_numer_denom()
    
    # We manually format the denominator for clarity, extracting the coefficients.
    # The denominator has the form C*r*(k1*sqrt(k2) + k3)
    # This ensures all numbers are explicitly printed as requested.
    den_factor = den / (C*r)
    k1 = den_factor.coeff(sympy.sqrt(3))
    k3 = den_factor.subs(sympy.sqrt(3), 0)
    
    print("The cutoff angular frequency w_c is:")
    print(f"w_c = {num} / (C * r * ({k1}*sqrt(3) + ({k3})))")


solve_cutoff_frequency()