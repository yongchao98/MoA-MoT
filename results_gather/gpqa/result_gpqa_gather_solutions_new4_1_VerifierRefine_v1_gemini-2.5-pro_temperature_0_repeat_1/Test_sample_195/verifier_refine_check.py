import sympy
from sympy import sqrt, Symbol, Eq, solve, limit, oo

def check_correctness():
    """
    Verifies the formula for the maximum speed of a relativistic harmonic oscillator.
    
    The function performs two checks:
    1. Derives the formula from the principle of conservation of energy and compares it to Option D.
    2. Checks if the formula from Option D correctly reduces to the classical formula (Option C)
       in the non-relativistic limit (c -> infinity).
    """
    try:
        # Define the physical symbols
        m = Symbol('m', positive=True, real=True)
        k = Symbol('k', positive=True, real=True)
        A = Symbol('A', positive=True, real=True)
        c = Symbol('c', positive=True, real=True)
        v_max = Symbol('v_max', real=True)

        # --- Check 1: Symbolic Derivation from First Principles ---

        # Total energy at maximum amplitude (x=A, v=0)
        # E_total = (rest energy) + (potential energy)
        energy_at_max_amplitude = m * c**2 + (k * A**2) / 2
        
        # Total energy at equilibrium (x=0, v=v_max)
        # E_total = (relativistic energy at v_max)
        gamma_max = 1 / sqrt(1 - v_max**2 / c**2)
        energy_at_equilibrium = gamma_max * m * c**2
        
        # Apply conservation of energy: E(x=A) = E(x=0)
        conservation_equation = Eq(energy_at_max_amplitude, energy_at_equilibrium)
        
        # Solve the equation for v_max
        solutions = solve(conservation_equation, v_max)
        
        # The solver returns two solutions (+v_max and -v_max). We take the positive one.
        derived_expression = None
        for sol in solutions:
            # The positive solution will not start with a minus sign when printed
            if not str(sol).strip().startswith('-'):
                derived_expression = sol
                break
        
        if derived_expression is None:
            return "Failed to isolate a positive solution for v_max from the conservation equation."

        # The formula from the proposed correct answer (Option D)
        option_D_expression = c * sqrt(1 - 1 / (1 + (k * A**2) / (2 * m * c**2))**2)

        # Check if the derived expression is equivalent to Option D
        # sympy.simplify(expr1 - expr2) == 0 is a robust way to check for equivalence
        if sympy.simplify(derived_expression - option_D_expression) != 0:
            return (f"The derived formula does not match Option D.\n"
                    f"Derived: {derived_expression}\n"
                    f"Option D: {option_D_expression}")

        # --- Check 2: Verification of the Classical Limit ---

        # The classical formula for v_max (Option C)
        classical_expression = sqrt(k * A**2 / m)
        
        # Calculate the limit of the relativistic formula (Option D) as c -> infinity
        classical_limit_of_D = limit(option_D_expression, c, oo)
        
        # Check if the limit matches the classical formula
        if sympy.simplify(classical_limit_of_D - classical_expression) != 0:
            return (f"The classical limit (c -> oo) of Option D does not match the classical formula (Option C).\n"
                    f"Limit of Option D: {classical_limit_of_D}\n"
                    f"Option C: {classical_expression}")

        # --- Final Check on Other Options ---
        # Option A/B would lead to v_max > c, which is physically impossible.
        # Option C is only the classical limit, not the general relativistic solution.
        
        return "Correct"

    except Exception as e:
        return f"An error occurred during the verification process: {e}"

# Run the check
result = check_correctness()
print(result)