import math

def check_answer():
    """
    Checks the correctness of the LLM's answer for the physics problem.
    """
    # Define the options as provided in the question
    options = {
        'A': {'fraction': 1/4, 'dependence': -3},
        'B': {'fraction': 1/4, 'dependence': -4},
        'C': {'fraction': 1/2, 'dependence': -4},
        'D': {'fraction': 3/4, 'dependence': -6}
    }

    # The final answer provided by the LLM
    llm_answer_choice = 'B'

    # --- Step 1: Calculate the correct fraction based on physics principles ---
    # The power radiated per unit solid angle (dP/dΩ) for an electric dipole
    # oscillating along the z-axis is proportional to sin²(θ).
    # The maximum power (A) is at θ = 90°, where sin²(90°) = 1.
    # We need to find the fraction at θ = 30°.
    theta_degrees = 30
    theta_radians = math.radians(theta_degrees)
    
    # The fraction is sin²(30°) / sin²(90°)
    correct_fraction = (math.sin(theta_radians))**2 / (math.sin(math.radians(90)))**2

    # --- Step 2: Determine the correct wavelength dependence ---
    # The power radiated by an oscillating electric dipole is proportional to ω⁴.
    # Since ω (angular frequency) is proportional to 1/λ (wavelength),
    # the power is proportional to (1/λ)⁴ = λ⁻⁴.
    correct_dependence_power = -4

    # --- Step 3: Check if the LLM's chosen option matches the calculated correct values ---
    if llm_answer_choice not in options:
        return f"Invalid option '{llm_answer_choice}'. The available options are A, B, C, D."

    chosen_option_data = options[llm_answer_choice]
    
    is_fraction_correct = math.isclose(chosen_option_data['fraction'], correct_fraction)
    is_dependence_correct = (chosen_option_data['dependence'] == correct_dependence_power)

    if is_fraction_correct and is_dependence_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_fraction_correct:
            error_messages.append(f"The fraction of power is incorrect. The correct fraction is sin²(30°) = {correct_fraction:.2f}, but option {llm_answer_choice} has a fraction of {chosen_option_data['fraction']:.2f}.")
        if not is_dependence_correct:
            error_messages.append(f"The wavelength dependence is incorrect. The correct dependence is λ^({correct_dependence_power}), but option {llm_answer_choice} has a dependence of λ^({chosen_option_data['dependence']}).")
        
        # Find the correct option for a more detailed error message
        correct_key = None
        for key, value in options.items():
            if math.isclose(value['fraction'], correct_fraction) and value['dependence'] == correct_dependence_power:
                correct_key = key
                break
        
        if correct_key:
            error_messages.append(f"The LLM chose option {llm_answer_choice}, but the correct option is {correct_key}.")
        else:
            error_messages.append("None of the provided options match the physically correct answer.")

        return "\n".join(error_messages)

# Run the check and print the result
result = check_answer()
print(result)