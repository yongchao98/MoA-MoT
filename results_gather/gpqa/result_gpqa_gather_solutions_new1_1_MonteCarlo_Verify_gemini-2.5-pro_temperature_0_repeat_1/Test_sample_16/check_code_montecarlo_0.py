import math

def check_chemistry_equilibrium():
    """
    This function verifies the answer to the Ca-EDTA equilibrium problem.
    It calculates the concentration of Ca2+ from first principles and compares
    it to the value given in the selected answer option.
    """
    # --- Define problem constants ---
    initial_complex_conc = 0.02  # M
    K_formation = 5e10

    # --- Define the multiple-choice options from the question ---
    options = {
        'A': 6.3e-7,
        'B': 2.0e-2,
        'C': 1.0e-2,
        'D': 5.0e-3
    }

    # --- The final answer provided by the LLM ---
    llm_selected_option = 'A'

    # --- Step 1: Calculate the dissociation constant (Kd) ---
    # The dissociation reaction is the reverse of the formation reaction.
    K_dissociation = 1 / K_formation

    # --- Step 2: Solve the equilibrium equation for [Ca2+] (x) ---
    # The equilibrium expression is: Kd = x^2 / (initial_conc - x)
    # This rearranges to a quadratic equation: x^2 + Kd*x - Kd*initial_conc = 0
    # We solve for x using the quadratic formula: x = [-b ± sqrt(b^2 - 4ac)] / 2a
    # where a=1, b=Kd, c=-Kd*initial_conc
    a = 1
    b = K_dissociation
    c = -K_dissociation * initial_complex_conc

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        return "Calculation Error: The discriminant is negative, so no real solution exists."

    # The concentration 'x' must be positive, so we take the positive root.
    calculated_ca_conc = (-b + math.sqrt(discriminant)) / (2 * a)

    # --- Step 3: Verify the LLM's answer ---
    # Check if the selected option is valid
    if llm_selected_option not in options:
        return f"Invalid Answer: The selected option '{llm_selected_option}' is not one of the possible choices {list(options.keys())}."

    # Get the value corresponding to the selected option
    answer_value = options[llm_selected_option]

    # Compare the calculated value with the answer's value using a relative tolerance.
    # A tolerance is needed because the option value is rounded.
    if math.isclose(calculated_ca_conc, answer_value, rel_tol=0.01): # 1% relative tolerance
        return "Correct"
    else:
        return (f"Incorrect. The calculated concentration of Ca2+ is approximately {calculated_ca_conc:.3e} M. "
                f"The selected answer '{llm_selected_option}' corresponds to a value of {answer_value:.3e} M. "
                f"The calculated value does not match the answer's value within a reasonable tolerance.")

# Run the check
result = check_chemistry_equilibrium()
print(result)