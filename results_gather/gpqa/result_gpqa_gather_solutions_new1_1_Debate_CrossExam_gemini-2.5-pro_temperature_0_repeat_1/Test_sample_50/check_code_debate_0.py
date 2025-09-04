import sympy

def check_potential_energy_answer():
    """
    Checks the correctness of the selected answer for the potential energy of a charge
    near a grounded conducting sphere.
    """
    # Define symbolic variables for the physical quantities.
    # k: Coulomb's constant, q: charge, R: sphere radius, d: distance from center.
    # All are assumed to be positive real numbers, and d > R.
    k, q, R, d = sympy.symbols('k q R d', real=True, positive=True)

    # The final answer provided by the LLM to be checked.
    llm_final_answer = 'D'

    # Define the formulas for each multiple-choice option as given in the question.
    options = {
        'A': -(1/2) * k * q**2 * R**2 / (d**2 - R**2),
        'B': -k * q**2 * d / (d**2 - R**2),
        'C': -(1/2) * k * q**2 * d / (d**2 + R**2),
        'D': -(1/2) * k * q**2 * R / (d**2 - R**2)
    }

    # --- Derivation of the correct formula using the method of images ---
    
    # 1. Define the image charge (q_prime) and its position (b).
    q_prime = -q * (R / d)
    b = R**2 / d

    # 2. Calculate the distance between the real charge q and the image charge q_prime.
    distance_q_q_prime = d - b

    # 3. Calculate the interaction energy between the real charge and the image charge.
    U_interaction = k * q * q_prime / distance_q_q_prime

    # 4. The net potential energy of the system is half of the interaction energy.
    U_correct_derived = (1/2) * U_interaction

    # Simplify the derived expression for a clean comparison.
    U_correct_simplified = sympy.simplify(U_correct_derived)

    # --- Verification ---

    # Get the formula corresponding to the LLM's final answer.
    llm_formula = options.get(llm_final_answer)
    if llm_formula is None:
        return f"Invalid answer choice '{llm_final_answer}'. The options are A, B, C, D."

    # Compare the LLM's chosen formula with the correctly derived formula.
    # sympy.simplify(expr1 - expr2) == 0 is a robust way to check for symbolic equality.
    if sympy.simplify(llm_formula - U_correct_simplified) == 0:
        return "Correct"
    else:
        # If the answer is wrong, find the correct option and provide a reason.
        correct_option = None
        for option, formula in options.items():
            if sympy.simplify(formula - U_correct_simplified) == 0:
                correct_option = option
                break
        
        reason = (f"The final answer is {llm_final_answer}, which corresponds to the formula: {llm_formula}.\n"
                  f"However, the correct formula derived from the method of images is: {U_correct_simplified}.\n")
        
        if correct_option:
            reason += f"The correct formula matches option {correct_option}."
        else:
            reason += "None of the given options match the correct formula."
            
        return reason

# Run the check and print the result.
result = check_potential_energy_answer()
print(result)