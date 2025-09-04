import sympy
import re

def check_final_answer_correctness():
    """
    Checks the correctness of the final LLM answer for the given physics problem.

    The function performs the following steps:
    1.  Defines the correct physical formula for the potential energy using the method of images.
    2.  Parses the multiple-choice options from the original question.
    3.  Parses the final selected answer (e.g., 'D') from the provided LLM response.
    4.  Compares the formula corresponding to the selected answer with the correct physical formula.
    5.  If the answer is correct, it returns "Correct".
    6.  If incorrect, it returns a detailed reason, including an analysis of physical constraints.
    """
    # Define symbolic variables. Using k_c for Coulomb's constant.
    # All variables are positive real numbers. We assume d > R for the charge to be outside the sphere.
    k_c, q, R, d = sympy.symbols('k_c q R d', real=True, positive=True)

    # --- Step 1: Define the correct physical formula from first principles ---
    # The net potential energy of the system is half the interaction energy
    # between the real charge q and its image charge q' = -qR/d.
    # U = (1/2) * U_interaction = (1/2) * [k * q * q' / (d - R^2/d)]
    correct_formula = -sympy.Rational(1, 2) * k_c * q**2 * R / (d**2 - R**2)

    # --- Step 2: Define the options from the original question ---
    # Note: The question uses 'k', which we map to our symbolic 'k_c'.
    options = {
        'A': -k_c * q**2 * d / (d**2 - R**2),
        'B': -sympy.Rational(1, 2) * k_c * q**2 * R**2 / (d**2 - R**2),
        'C': -sympy.Rational(1, 2) * k_c * q**2 * d / (d**2 + R**2),
        'D': -sympy.Rational(1, 2) * k_c * q**2 * R / (d**2 - R**2)
    }

    # --- Step 3: Identify the final answer from the provided text ---
    # The final answer to be checked is the one from the last "meta-analysis" block.
    final_answer_text = """
    Based on an analysis of the candidate answers and a first-principles derivation, here is a step-by-step breakdown of the problem.
    ...
    <<<D>>>
    """ # Abridged for clarity
    
    # The final answer is D.
    chosen_answer_key = 'D'
    chosen_formula = options.get(chosen_answer_key)

    # --- Step 4: Compare the chosen answer with the correct formula and constraints ---
    if sympy.simplify(chosen_formula - correct_formula) == 0:
        # The formula is symbolically correct. As a final check, ensure it also passes physical constraints.
        # Constraint 1: Limit as d -> R+ should be -infinity.
        limit_at_R = sympy.limit(chosen_formula, d, R, dir='+')
        # Constraint 2: Limit as d -> infinity should fall off as 1/d^2.
        # We check this by seeing if limit(formula * d^2) is a finite constant.
        limit_at_inf = sympy.limit(chosen_formula * d**2, d, sympy.oo)

        if limit_at_R == -sympy.oo and limit_at_inf.is_finite and limit_at_inf != 0:
            return "Correct"
        else:
            # This case should not be reached if the formula is symbolically correct, but it's a good safeguard.
            return f"Incorrect. The formula for {chosen_answer_key} is symbolically correct but fails physical constraint checks."
    else:
        # The chosen answer is incorrect. Provide a detailed reason.
        reason = f"Incorrect. The chosen answer is {chosen_answer_key}, which corresponds to the formula: {chosen_formula}.\n"
        reason += f"The correct formula, derived from the method of images, is: {correct_formula}.\n\n"
        reason += f"Analysis of the chosen formula for answer {chosen_answer_key}:\n"

        # Constraint Check 1: Factor of 1/2
        if 'Rational(1, 2)' not in str(sympy.srepr(chosen_formula)):
            reason += "- FAILED: The formula is missing the crucial factor of 1/2.\n"
        
        # Constraint Check 2: Numerator (proportional to R, not d or R^2)
        if d in chosen_formula.as_numer_denom()[0].free_symbols:
             reason += "- FAILED: The numerator should be proportional to R, not d.\n"
        elif R**2 in chosen_formula.as_numer_denom()[0].free_symbols:
             reason += "- FAILED: The numerator should be proportional to R, not R^2 (incorrect dimensions).\n"

        # Constraint Check 3: Denominator (d^2 - R^2)
        if (d**2 + R**2) in chosen_formula.as_numer_denom()[1].args:
            limit_at_R = sympy.limit(chosen_formula, d, R, dir='+')
            reason += f"- FAILED: The denominator is (d^2 + R^2). As d->R, the energy should go to -infinity, but it approaches a finite value of {limit_at_R}.\n"

        # Constraint Check 4: Asymptotic behavior (should be ~1/d^2)
        order_at_inf = sympy.O(chosen_formula, (d, sympy.oo))
        if order_at_inf != sympy.O(1/d**2):
            reason += f"- FAILED: As d->infinity, the energy should fall off as 1/d^2. The chosen formula's behavior is {order_at_inf}.\n"

        return reason.strip()

# Run the check
result = check_final_answer_correctness()
print(result)