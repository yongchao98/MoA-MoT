import math

def check_physics_answer():
    """
    Checks the correctness of the answer based on the physics of electric dipole radiation.
    """
    # Define the options from the question
    # Format: { 'Option_Letter': (Fraction, Power_of_Lambda) }
    options = {
        'A': (3/4, -6),
        'B': (1/4, -3),
        'C': (1/4, -4),
        'D': (1/2, -4)
    }

    # The final answer provided by the LLM to be checked
    llm_answer_key = 'C'

    # --- Step 1: Determine the correct values based on physics ---

    # Constraint 1: Angular Dependence (Fraction of Maximum Power)
    # For an electric dipole, power is proportional to sin^2(θ).
    # Maximum power is at θ = 90 degrees. We need the fraction at θ = 30 degrees.
    # Fraction = sin^2(30°) / sin^2(90°) = (0.5)^2 / 1^2 = 0.25
    angle_deg = 30
    correct_fraction = math.sin(math.radians(angle_deg))**2

    # Constraint 2: Wavelength Dependence
    # Radiated power P is proportional to ω^4 (angular frequency).
    # Since ω ∝ 1/λ (wavelength), P ∝ (1/λ)^4 = λ^(-4).
    correct_lambda_power = -4

    # --- Step 2: Check if the LLM's chosen answer satisfies the constraints ---

    if llm_answer_key not in options:
        return f"Invalid answer key '{llm_answer_key}'. The key must be one of {list(options.keys())}."

    llm_fraction, llm_lambda_power = options[llm_answer_key]

    # Use math.isclose for robust floating-point comparison
    is_fraction_correct = math.isclose(llm_fraction, correct_fraction)
    is_lambda_power_correct = (llm_lambda_power == correct_lambda_power)

    if is_fraction_correct and is_lambda_power_correct:
        return "Correct"
    else:
        # If the answer is wrong, explain why.
        error_messages = []
        if not is_fraction_correct:
            error_messages.append(
                f"The fraction of maximum power is incorrect. "
                f"Based on sin^2(30°), the expected fraction is {correct_fraction}, "
                f"but option {llm_answer_key} provides {llm_fraction}."
            )
        if not is_lambda_power_correct:
            error_messages.append(
                f"The wavelength dependence is incorrect. "
                f"Based on dipole radiation theory, the power should be proportional to λ^({correct_lambda_power}), "
                f"but option {llm_answer_key} provides λ^({llm_lambda_power})."
            )
        
        # Find the truly correct option for a more complete error message
        truly_correct_key = None
        for key, (frac, power) in options.items():
            if math.isclose(frac, correct_fraction) and power == correct_lambda_power:
                truly_correct_key = key
                break
        
        if truly_correct_key and truly_correct_key != llm_answer_key:
             error_messages.append(f"The correct option is {truly_correct_key}.")

        return "Incorrect. " + " ".join(error_messages)

# Run the check and print the result
result = check_physics_answer()
print(result)