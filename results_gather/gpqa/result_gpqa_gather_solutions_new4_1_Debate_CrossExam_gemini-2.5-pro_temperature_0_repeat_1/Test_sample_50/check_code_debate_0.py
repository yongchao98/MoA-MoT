import sympy

def check_potential_energy_answer():
    """
    Checks the correctness of the final answer for the potential energy of a charge
    near a grounded conducting sphere.

    The function symbolically derives the correct formula using the method of images
    and compares it with the formula corresponding to the given answer 'A'.
    """
    try:
        # Define symbolic variables for the physical quantities.
        # 'k' is Coulomb's constant.
        # 'q' is the charge.
        # 'R' is the radius of the sphere.
        # 'd' is the distance of the charge from the center of the sphere.
        # All are assumed to be positive real numbers.
        k, q, R, d = sympy.symbols('k q R d', real=True, positive=True)

        # The physical constraint is that the charge is outside the sphere, so d > R.
        # This ensures the denominator (d**2 - R**2) is positive and the formula is well-defined.

        # --- Step 1: Derive the correct formula from first principles (method of images) ---
        
        # Image charge magnitude and position
        q_prime = -q * (R / d)
        b = R**2 / d
        
        # Distance between the real charge q and the image charge q'
        distance_q_q_prime = d - b
        
        # Potential at the location of q due to the image charge q'
        V_induced = k * q_prime / distance_q_q_prime
        
        # The potential energy of the system is (1/2) * q * V_induced
        correct_formula = sympy.Rational(1, 2) * q * V_induced
        
        # Simplify the derived formula
        correct_formula_simplified = sympy.simplify(correct_formula)

        # --- Step 2: Define the formulas from the multiple-choice options ---
        
        # A) U = - (1/2) * k * q^2 * R / (d^2 - R^2)
        option_A = - (sympy.Rational(1, 2) * k * q**2 * R) / (d**2 - R**2)
        
        # B) U = - k * q^2 * d / (d^2 - R^2)
        option_B = - (k * q**2 * d) / (d**2 - R**2)
        
        # C) U = - (1/2) * k * q^2 * R^2 / (d^2 - R^2)
        option_C = - (sympy.Rational(1, 2) * k * q**2 * R**2) / (d**2 - R**2)
        
        # D) U = - (1/2) * k * q^2 * d / (d^2 + R^2)
        option_D = - (sympy.Rational(1, 2) * k * q**2 * d) / (d**2 + R**2)

        # --- Step 3: Compare the given answer with the correct formula ---
        
        # The final answer provided by the LLM is 'A'.
        llm_selected_formula = option_A

        # Check if the LLM's chosen formula is symbolically equivalent to the correct one.
        # sympy.simplify(expr1 - expr2) == 0 is a robust way to check for equality.
        if sympy.simplify(llm_selected_formula - correct_formula_simplified) == 0:
            # The analysis in the prompt is also sound. It correctly identifies the
            # method of images, derives the correct formula, and matches it to option A.
            return "Correct"
        else:
            reason = "The final answer 'A' is incorrect.\n"
            reason += f"The formula for option A is: {llm_selected_formula}\n"
            reason += f"The correct formula derived from the method of images is: {correct_formula_simplified}\n"
            reason += "The selected formula does not match the correct physical derivation."
            return reason
            
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result.
result = check_potential_energy_answer()
print(result)