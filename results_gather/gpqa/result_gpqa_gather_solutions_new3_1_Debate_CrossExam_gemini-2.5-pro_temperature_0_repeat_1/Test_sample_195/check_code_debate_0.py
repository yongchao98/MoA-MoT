import sympy

def check_answer():
    """
    Checks the correctness of the answer for the relativistic harmonic oscillator problem.
    It uses symbolic mathematics to derive the solution from first principles and
    compares it to the provided answer (Option C).
    """
    try:
        # 1. Define all symbolic variables.
        # We assume all physical quantities are positive real numbers.
        m, k, A, c = sympy.symbols('m k A c', positive=True, real=True)
        v_max = sympy.symbols('v_max', positive=True, real=True)

        # The final answer provided in the prompt is C.
        # Let's define the expression for option C to check against.
        # C) v_max = c * sqrt(1 - 1 / (1 + kA^2/(2mc^2))^2)
        option_C_expr = c * sympy.sqrt(1 - 1 / (1 + k * A**2 / (2 * m * c**2))**2)

        # 2. Set up the conservation of energy equation from first principles.
        # Total energy E_total = gamma * m * c^2 + (1/2) * k * x^2

        # State 1: Energy at maximum amplitude (x=A, v=0)
        # At v=0, the Lorentz factor gamma is 1.
        E_at_amplitude = m * c**2 + sympy.Rational(1, 2) * k * A**2

        # State 2: Energy at equilibrium (x=0, v=v_max)
        # At x=0, the potential energy is 0.
        gamma_max = 1 / sympy.sqrt(1 - v_max**2 / c**2)
        E_at_equilibrium = gamma_max * m * c**2

        # 3. Equate the energies to represent the conservation principle.
        energy_conservation_eq = sympy.Eq(E_at_amplitude, E_at_equilibrium)

        # 4. Solve the equation for v_max.
        # sympy.solve returns a list of solutions. We need the positive one.
        solutions = sympy.solve(energy_conservation_eq, v_max)
        
        derived_solution = None
        for sol in solutions:
            # The speed v_max must be a positive quantity.
            # We can check if the expression is guaranteed to be positive given our assumptions.
            if sol.is_positive:
                derived_solution = sol
                break
        
        if derived_solution is None:
            return "Error: Could not derive a valid positive solution for v_max from the conservation of energy principle."

        # 5. Compare the derived solution with the expression from option C.
        # sympy.simplify(expr1 - expr2) will be 0 if they are algebraically identical.
        if sympy.simplify(derived_solution - option_C_expr) == 0:
            # The derived solution matches option C. This confirms the answer is correct.
            # As an additional check, we can verify the non-relativistic (classical) limit.
            
            # 6. Sanity Check: Classical Limit (c -> infinity)
            # The expected classical result is v_max = sqrt(k*A^2/m)
            classical_v_max_expr = sympy.sqrt(k * A**2 / m)
            
            # Compute the limit of our derived solution as c approaches infinity.
            limit_of_derived_solution = sympy.limit(derived_solution, c, sympy.oo)
            
            if sympy.simplify(limit_of_derived_solution - classical_v_max_expr) == 0:
                # The classical limit is also correct.
                return "Correct"
            else:
                return f"Incorrect. The provided answer C is algebraically correct, but it fails the classical limit check. The limit of the derived expression is {limit_of_derived_solution}, while the expected classical result is {classical_v_max_expr}."
        else:
            return f"Incorrect. The provided answer is C, but the correct derivation from the principle of conservation of relativistic energy leads to a different expression: {derived_solution}."

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check and print the result.
result = check_answer()
print(result)