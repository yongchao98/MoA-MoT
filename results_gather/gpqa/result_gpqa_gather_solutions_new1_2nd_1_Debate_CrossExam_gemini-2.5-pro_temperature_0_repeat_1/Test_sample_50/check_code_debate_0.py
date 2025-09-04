import sympy

def check_correctness():
    """
    Checks the correctness of the LLM's answer by deriving the correct formula
    and comparing it to the chosen option. It also verifies the answer against
    key physical constraints.
    """
    # Define symbolic variables for the physical quantities.
    # 'd > R' is a necessary physical constraint for the charge to be outside the sphere.
    k, q, R, d = sympy.symbols('k q R d', positive=True, real=True)

    # The final answer provided by the LLM.
    llm_answer_str = "<<<A>>>"
    
    # Extract the letter from the answer string.
    try:
        selected_option = llm_answer_str.strip().replace('<', '').replace('>', '')
    except:
        return "Invalid answer format. Could not extract option letter."

    # Define the formulas for each option as given in the question.
    options = {
        'A': - (sympy.S(1)/2) * k * q**2 * R / (d**2 - R**2),
        'B': - (sympy.S(1)/2) * k * q**2 * d / (d**2 + R**2),
        'C': - (sympy.S(1)/2) * k * q**2 * R**2 / (d**2 - R**2),
        'D': - k * q**2 * d / (d**2 - R**2)
    }

    # --- Step 1: Derive the correct formula from first principles (Method of Images) ---
    # Image charge q' = -q*(R/d) at position b = R**2/d
    # Interaction energy U_int = k*q*q' / (d-b) = -k*q**2*R / (d**2 - R**2)
    # The net potential energy of the system is half the interaction energy.
    correct_formula = (sympy.S(1)/2) * (-k * q**2 * R / (d**2 - R**2))
    
    # --- Step 2: Check if the selected option is valid ---
    if selected_option not in options:
        return f"Invalid option '{selected_option}'. Please choose from A, B, C, D."

    llm_formula = options[selected_option]

    # --- Step 3: Compare the LLM's formula with the correct one ---
    # sympy.simplify will return 0 if the expressions are equivalent.
    if sympy.simplify(llm_formula - correct_formula) == 0:
        return "Correct"
    else:
        # --- Step 4: If incorrect, provide a reason based on physical constraints ---
        
        # Constraint 1: Asymptotic behavior (as d -> infinity, U should fall as 1/d^2)
        # We check the leading term of the series expansion for large d.
        series_expansion = sympy.series(llm_formula, d, sympy.oo, n=3).removeO()
        if not series_expansion.as_leading_term(d).has(d**-2):
            return (f"Incorrect. The formula for option {selected_option} fails the asymptotic behavior test. "
                    f"As d -> infinity, the potential energy should fall off as 1/d^2. "
                    f"The provided formula {llm_formula} falls off as {series_expansion.as_leading_term(d)}.")

        # Constraint 2: Boundary condition (as d -> R+, U should go to -infinity)
        # This requires the denominator to go to zero.
        denominator = sympy.fraction(llm_formula)[1]
        if sympy.limit(denominator, d, R, dir='+') != 0:
            return (f"Incorrect. The formula for option {selected_option} fails the boundary condition test. "
                    f"As the charge approaches the sphere (d -> R+), the potential energy should go to -infinity. "
                    f"The provided formula {llm_formula} approaches a finite value.")

        # Constraint 3: Check for the crucial factor of 1/2.
        # The interaction energy is twice the system energy.
        interaction_energy = -k * q**2 * R / (d**2 - R**2)
        if sympy.simplify(llm_formula - interaction_energy) == 0:
            return (f"Incorrect. The formula for option {selected_option} is missing the crucial factor of 1/2. "
                    f"It represents the interaction energy, not the net potential energy of the system.")
        
        # General failure message if the formula is wrong in another way (e.g., wrong power of R).
        return f"Incorrect. The formula for option {selected_option}, which is {llm_formula}, does not match the correct physical formula {correct_formula}."

# Run the check
result = check_correctness()
print(result)