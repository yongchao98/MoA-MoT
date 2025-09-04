import math

def check_answer():
    """
    Checks if the selected answer's formula matches the correct physical formula.
    """
    # The final answer provided by the analysis is 'D'.
    final_answer_option = 'D'

    # Define a set of sample physical values.
    # We must satisfy the constraint d > R.
    k = 8.98755e9  # Coulomb's constant
    q = 2e-9       # Charge in Coulombs
    R = 0.1        # Radius in meters
    d = 0.2        # Distance in meters

    # --- Formulas Definition ---

    # This is the known correct formula from physics textbooks.
    def get_correct_potential_energy(k_val, q_val, R_val, d_val):
        """Calculates U = - (1/2) * k * q^2 * R / (d^2 - R^2)"""
        if d_val <= R_val:
            raise ValueError("Distance 'd' must be greater than radius 'R'.")
        numerator = -0.5 * k_val * q_val**2 * R_val
        denominator = d_val**2 - R_val**2
        return numerator / denominator

    # Define a dictionary mapping options to their respective formulas.
    formulas = {
        'A': lambda k, q, R, d: -k * q**2 * d / (d**2 - R**2),
        'B': lambda k, q, R, d: -0.5 * k * q**2 * R**2 / (d**2 - R**2),
        'C': lambda k, q, R, d: -0.5 * k * q**2 * d / (d**2 + R**2),
        'D': lambda k, q, R, d: -0.5 * k * q**2 * R / (d**2 - R**2)
    }

    # --- Verification Logic ---

    # Get the function corresponding to the final answer.
    chosen_formula = formulas.get(final_answer_option)

    if not chosen_formula:
        return f"Invalid option '{final_answer_option}' selected. Cannot perform check."

    # Calculate the energy using the correct formula and the chosen formula.
    try:
        correct_value = get_correct_potential_energy(k, q, R, d)
        chosen_value = chosen_formula(k, q, R, d)
    except (ValueError, ZeroDivisionError) as e:
        return f"Error during calculation: {e}"

    # Compare the results. Use math.isclose for robust floating-point comparison.
    if math.isclose(correct_value, chosen_value):
        return "Correct"
    else:
        # If they don't match, find out which option *is* correct.
        correct_option = None
        for option, func in formulas.items():
            if math.isclose(func(k, q, R, d), correct_value):
                correct_option = option
                break
        
        reason = (f"Incorrect. The final answer claims option '{final_answer_option}' is correct, "
                  f"but its formula is wrong.\n"
                  f"The correct formula is U = - (1/2) * k * q^2 * R / (d^2 - R^2), "
                  f"which corresponds to option '{correct_option}'.")
        return reason

# Run the check and print the result.
result = check_answer()
print(result)