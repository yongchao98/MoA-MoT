import sympy

def check_relativistic_oscillator_answer():
    """
    Verifies the correct formula for the maximum speed of a relativistic harmonic oscillator.

    The function performs three checks:
    1. Derives the correct formula from first principles using symbolic math.
    2. Checks if the derived formula matches the proposed answer (Option C).
    3. Performs a sanity check by taking the non-relativistic limit (c -> infinity)
       and ensuring it matches the classical result (Option A).
    4. Performs a physical plausibility check to ensure v_max is always less than c.
    """
    # 1. Define symbolic variables
    # All physical quantities are positive real numbers
    m, k, A, c = sympy.symbols('m k A c', positive=True, real=True)
    v_max = sympy.symbols('v_max', real=True)

    # 2. Derive the solution from the principle of conservation of energy
    
    # Total relativistic energy E_total = gamma*m*c**2 + U(x)
    # At maximum amplitude (x=A), velocity is 0, so gamma=1.
    E_at_amplitude = m * c**2 + sympy.Rational(1, 2) * k * A**2

    # At equilibrium (x=0), potential energy is 0, and speed is v_max.
    gamma_max = 1 / sympy.sqrt(1 - v_max**2 / c**2)
    E_at_equilibrium = gamma_max * m * c**2

    # Equate the energies to represent conservation
    energy_conservation_eq = sympy.Eq(E_at_amplitude, E_at_equilibrium)

    # Solve for v_max. It's easier to solve for gamma_max first, then v_max.
    gamma_max_expr = sympy.solve(energy_conservation_eq, gamma_max)[0]
    
    # Now solve for v_max from the definition of gamma_max
    # gamma_max = 1 / sqrt(1 - v_max**2 / c**2)
    v_max_solutions = sympy.solve(sympy.Eq(gamma_max_expr, 1 / sympy.sqrt(1 - v_max**2 / c**2)), v_max)
    
    # The solutions will be positive and negative; speed is the positive value.
    derived_solution = next(sol for sol in v_max_solutions if sol.is_positive)

    # 3. Define the candidate answers symbolically
    option_A = sympy.sqrt(k * A**2 / m)
    option_C = c * sympy.sqrt(1 - 1 / (1 + (k * A**2) / (2 * m * c**2))**2)
    
    # 4. Check if the derived solution matches the proposed answer (C)
    # sympy.simplify will return 0 if the expressions are algebraically equivalent.
    if sympy.simplify(derived_solution - option_C) != 0:
        return (f"Incorrect. The derivation from first principles yields:\n{derived_solution}\n"
                f"This is not algebraically equivalent to the proposed answer C:\n{option_C}")

    # 5. Sanity Check 1: Non-relativistic limit (c -> infinity)
    # The result should be the classical formula (Option A)
    classical_limit = sympy.limit(option_C, c, sympy.oo)
    if sympy.simplify(classical_limit - option_A) != 0:
        return (f"Incorrect. The non-relativistic limit (c -> oo) of option C is:\n{classical_limit}\n"
                f"This does not match the expected classical result (Option A):\n{option_A}")

    # 6. Sanity Check 2: Physical Plausibility (v_max < c)
    # For any positive real values of m, k, A, c, v_max must be less than c.
    # This is equivalent to checking if option_C / c < 1.
    # The term X = k*A**2/(2*m*c**2) is > 0.
    # The term 1/(1+X)**2 is between 0 and 1.
    # The term inside the sqrt, 1 - (term between 0 and 1), is also between 0 and 1.
    # So the sqrt is between 0 and 1, which means v_max < c.
    # We can ask sympy to prove this inequality.
    is_slower_than_light = sympy.simplify(option_C < c)
    if is_slower_than_light is not sympy.true:
        return (f"Incorrect. The formula for option C does not satisfy the physical constraint v_max < c. "
                f"The check `option_C < c` evaluates to {is_slower_than_light} instead of True.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_relativistic_oscillator_answer()
print(result)