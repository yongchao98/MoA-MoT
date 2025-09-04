import sympy

def check_potential_energy_answer():
    """
    Checks the correctness of the LLM's answer for the potential energy of a charge near a grounded sphere.
    """
    # Define symbolic variables
    k, q, R, d = sympy.symbols('k q R d', real=True, positive=True)
    
    # The correct formula for the potential energy derived from the method of images.
    # The energy of the system is half the interaction energy between the real charge and the image charge.
    # Image charge q' = -q*R/d at position b = R^2/d.
    # Interaction energy U_int = k*q*q' / (d-b) = k*q*(-q*R/d) / (d - R**2/d) = -k*q**2*R / (d**2 - R**2)
    # System energy U = 1/2 * U_int
    correct_formula = -sympy.S(1)/2 * k * q**2 * R / (d**2 - R**2)

    # The options provided in the question
    options = {
        'A': -k * q**2 * d / (d**2 - R**2),
        'B': -sympy.S(1)/2 * k * q**2 * d / (d**2 + R**2),
        'C': -sympy.S(1)/2 * k * q**2 * R**2 / (d**2 - R**2),
        'D': -sympy.S(1)/2 * k * q**2 * R / (d**2 - R**2)
    }

    # The final answer provided by the LLM
    llm_answer_key = 'D'
    llm_answer_formula = options.get(llm_answer_key)

    if llm_answer_formula is None:
        return f"Invalid answer key '{llm_answer_key}' provided."

    # Check if the LLM's chosen formula is symbolically equal to the correct formula
    # sympy.simplify() is used to ensure the comparison is robust
    if sympy.simplify(correct_formula - llm_answer_formula) == 0:
        return "Correct"
    else:
        # If incorrect, explain the discrepancy
        reason = f"The provided answer corresponds to the formula: U = {llm_answer_formula}.\n"
        reason += f"The correct formula, derived from the method of images, is: U = {correct_formula}.\n"
        reason += "The answer is incorrect because it does not match the physically derived result."
        return reason

# Run the check
result = check_potential_energy_answer()
print(result)