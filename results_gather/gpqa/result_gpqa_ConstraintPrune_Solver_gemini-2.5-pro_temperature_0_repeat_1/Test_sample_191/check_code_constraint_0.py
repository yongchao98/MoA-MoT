import sympy

def check_conductor_field_answer():
    """
    Checks the correctness of the answer for the electric field outside a conductor with a cavity.

    The function verifies the answer based on two fundamental principles of electrostatics:
    1.  Electrostatic Shielding: The electric field outside a closed conductor is independent
        of the configuration of charges inside it. The formula must not depend on parameters
        describing the internal geometry (l, s, theta).
    2.  Shell Theorem (via Gauss's Law): The charge +q induced on the outer spherical surface
        distributes uniformly. The field at an external point is therefore identical to that
        of a point charge +q located at the center of the conductor. The formula must be
        proportional to q/L^2.
    """
    # The answer provided by the LLM
    llm_answer = 'D'

    # Define symbolic variables to represent the physical quantities
    # k represents the Coulomb constant 1/(4*pi*epsilon_o)
    k, q, R, s, r, L, l, theta = sympy.symbols('k q R s r L l theta', real=True, positive=True)

    # Define the options as symbolic expressions
    options = {
        'A': k * q / (l - s * sympy.cos(theta))**2,
        'B': k * q / (l + s * sympy.cos(theta))**2,
        'C': k * q / l**2,
        'D': k * q / L**2
    }

    # --- Constraint 1: Electrostatic Shielding ---
    # The correct formula must not depend on internal geometry variables: l, s, theta.
    # It should only depend on the distance from the conductor's center (L) and the charge (q).
    
    candidates_after_shielding = []
    for option_key, formula in options.items():
        variables_in_formula = formula.free_symbols
        # Check if any internal geometry variables are present
        if not any(var in variables_in_formula for var in [l, s, theta]):
            candidates_after_shielding.append(option_key)

    if llm_answer not in candidates_after_shielding:
        return (f"Incorrect. The answer '{llm_answer}' violates the principle of electrostatic shielding. "
                f"The electric field outside the conductor cannot depend on the internal parameters "
                f"l, s, or theta, which describe the position of the charge inside the cavity. "
                f"The conductor shields the exterior from the internal field configuration.")

    # --- Constraint 2: Shell Theorem / Gauss's Law ---
    # The induced charge +q on the outer surface is distributed uniformly.
    # The field outside a uniformly charged sphere is k*q / (distance_from_center)^2.
    # Here, the distance from the center is L.
    
    expected_formula = k * q / L**2
    correct_option = None

    for option_key in candidates_after_shielding:
        formula = options[option_key]
        # Check if the formula matches the expected form
        if sympy.simplify(formula - expected_formula) == 0:
            correct_option = option_key
            break # Found the correct option

    if correct_option is None:
        # This case is unlikely given the options, but included for completeness
        return "Incorrect. No option satisfies both electrostatic shielding and the shell theorem."

    if llm_answer == correct_option:
        return "Correct"
    else:
        # This case would trigger if the LLM passed shielding but failed the shell theorem check
        return (f"Incorrect. The answer '{llm_answer}' is independent of internal geometry, but its form "
                f"is incorrect. The correct formula, based on the shell theorem for the induced charge +q "
                f"on the outer surface, is proportional to q/L^2, which corresponds to option '{correct_option}'.")

# Run the check
result = check_conductor_field_answer()
print(result)