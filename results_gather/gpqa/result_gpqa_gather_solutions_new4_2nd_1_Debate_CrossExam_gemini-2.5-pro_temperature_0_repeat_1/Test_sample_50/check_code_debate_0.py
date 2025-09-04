import sympy

def check_answer_correctness():
    """
    Checks if the provided LLM answer for the potential energy of a charge
    near a grounded sphere is correct.
    """
    # Define symbolic variables for the physical quantities.
    # Using real=True and positive=True helps sympy with simplifications.
    k, q, R, d = sympy.symbols('k q R d', real=True, positive=True)

    # --- Ground Truth Formula Derivation ---
    # This is the standard result from the method of images.
    # 1. Image charge q_prime = -q * (R/d)
    # 2. Image position b = R**2 / d
    # 3. Distance between real and image charge = d - b = (d**2 - R**2) / d
    # 4. Potential at q due to image charge (V_induced) = k * q_prime / (d - b)
    #    V_induced = k * (-q * R / d) / ((d**2 - R**2) / d)
    #    V_induced = -k * q * R / (d**2 - R**2)
    # 5. Potential Energy U = (1/2) * q * V_induced
    ground_truth_formula = -(sympy.S(1)/2) * k * q**2 * R / (d**2 - R**2)

    # --- Formulas from the Multiple Choice Options ---
    # These are taken from the final analysis section of the provided answer.
    option_A = -(sympy.S(1)/2) * k * q**2 * d / (d**2 + R**2)
    option_B = -k * q**2 * d / (d**2 - R**2)
    option_C = -(sympy.S(1)/2) * k * q**2 * R / (d**2 - R**2)
    option_D = -(sympy.S(1)/2) * k * q**2 * R**2 / (d**2 - R**2)

    # The LLM's final answer is <<<C>>>
    llm_answer_letter = 'C'
    
    options_map = {
        'A': option_A,
        'B': option_B,
        'C': option_C,
        'D': option_D
    }
    
    llm_selected_formula = options_map.get(llm_answer_letter)

    # --- Verification Step ---
    # We check if the LLM's selected formula is symbolically equal to the ground truth.
    # A robust way to check for symbolic equality is to see if their difference simplifies to zero.
    if sympy.simplify(llm_selected_formula - ground_truth_formula) == 0:
        # The main physical constraint is that the charge is outside the sphere (d > R).
        # The formula U = - (1/2) * k * q^2 * R / (d^2 - R^2) is valid under this constraint.
        # Since d > R, d^2 - R^2 > 0, and the potential energy U is negative,
        # which is physically correct for an attractive force. The formula is consistent.
        return "Correct"
    else:
        return (f"Incorrect. The selected answer {llm_answer_letter} corresponds to the formula U = {llm_selected_formula}. "
                f"The correct formula, derived using the method of images, is U = {ground_truth_formula}.")

# Run the check and print the result.
result = check_answer_correctness()
print(result)