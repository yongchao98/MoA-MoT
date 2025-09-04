import sympy

def check_correctness():
    """
    Checks the correctness of the answer for the relativistic harmonic oscillator problem.
    
    The function performs three main checks on the proposed correct answer (Option C):
    1. Derives the solution from the principle of conservation of total relativistic energy and
       compares it to Option C.
    2. Checks if Option C correctly reduces to the classical result in the non-relativistic limit (c -> infinity).
    3. Verifies that Option C always satisfies the physical constraint v_max < c.
    """
    try:
        # 1. Define symbolic variables for the physical quantities.
        # All are defined as positive real numbers.
        m, k, A, c = sympy.symbols('m k A c', positive=True)
        v_max = sympy.Symbol('v_max', real=True)

        # 2. State the principle of conservation of total relativistic energy.
        # E_total = gamma * m * c**2 + 0.5 * k * x**2
        
        # Energy at maximum amplitude (x=A, v=0, gamma=1)
        E_at_amplitude = m * c**2 + sympy.Rational(1, 2) * k * A**2

        # Energy at equilibrium (x=0, v=v_max)
        gamma_max = 1 / sympy.sqrt(1 - v_max**2 / c**2)
        E_at_equilibrium = gamma_max * m * c**2

        # Formulate the conservation equation
        conservation_eq = sympy.Eq(E_at_amplitude, E_at_equilibrium)

        # 3. Solve the equation for v_max to get the derived solution.
        solutions = sympy.solve(conservation_eq, v_max)
        
        # The solver returns two solutions (+v_max and -v_max). We take the positive one for speed.
        derived_solution = None
        for sol in solutions:
            if sol.is_positive:
                derived_solution = sol
                break
        
        if derived_solution is None:
            return "Error: Could not derive a positive solution for v_max from the conservation of energy principle."

        # 4. Define the proposed correct answer (Option C) as a symbolic expression.
        option_C = c * sympy.sqrt(1 - 1 / (1 + k*A**2 / (2*m*c**2))**2)

        # 5. Check 1: Does the derived solution match Option C?
        if sympy.simplify(derived_solution - option_C) != 0:
            return (f"Incorrect: The proposed answer (Option C) does not match the result derived from first principles.\n"
                    f"Derived expression: {derived_solution}\n"
                    f"Option C: {option_C}")

        # 6. Check 2: Does Option C have the correct non-relativistic limit?
        # As c -> infinity, v_max should approach the classical result sqrt(k*A**2/m).
        classical_vmax = sympy.sqrt(k * A**2 / m)
        limit_of_C = sympy.limit(option_C, c, sympy.oo)

        if sympy.simplify(limit_of_C - classical_vmax) != 0:
            return (f"Incorrect: The non-relativistic limit (c -> oo) of Option C is incorrect.\n"
                    f"Expected limit: {classical_vmax}\n"
                    f"Actual limit: {limit_of_C}")

        # 7. Check 3: Does Option C satisfy the physical constraint v_max < c?
        # This is equivalent to checking if c**2 - v_max**2 is always positive.
        v_max_sq = option_C**2
        diff_sq = sympy.simplify(c**2 - v_max_sq)
        
        # .is_positive is a powerful tool in sympy to check the sign of an expression
        if not diff_sq.is_positive:
            return (f"Incorrect: Option C does not guarantee that v_max < c.\n"
                    f"The expression for c^2 - v_max^2 simplifies to {diff_sq}, which is not guaranteed to be positive.")

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_correctness()
print(result)