import math

def check_dipole_radiation_answer():
    """
    Checks the correctness of the answer based on the physics of electric dipole radiation.
    """
    # Define the options provided in the question
    # Format: { 'OptionLetter': (Fraction, PowerOfLambda) }
    options = {
        'A': (3/4, -6),
        'B': (1/4, -3),
        'C': (1/2, -4),
        'D': (1/4, -4)
    }

    # The final answer to be checked
    final_answer_letter = 'D'

    # --- Constraint 1: Calculate the expected fraction of power ---
    # The radiated power is proportional to sin²(θ).
    # The maximum is at θ=90°. The fraction at θ=30° is sin²(30°).
    theta_deg = 30
    expected_fraction = math.sin(math.radians(theta_deg))**2

    # --- Constraint 2: Determine the expected wavelength dependence ---
    # Radiated power is proportional to ω⁴, and ω is proportional to 1/λ.
    # Therefore, power is proportional to (1/λ)⁴ = λ⁻⁴.
    expected_power_of_lambda = -4

    # Retrieve the values from the chosen answer
    if final_answer_letter not in options:
        return f"Invalid option '{final_answer_letter}'. The options are A, B, C, D."

    chosen_fraction, chosen_power = options[final_answer_letter]

    # --- Verify the constraints ---
    is_fraction_correct = math.isclose(chosen_fraction, expected_fraction)
    is_power_correct = (chosen_power == expected_power_of_lambda)

    if is_fraction_correct and is_power_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_fraction_correct:
            error_messages.append(
                f"Constraint violated: Angular dependence. "
                f"The fraction of maximum power at 30° should be sin²(30°) = {expected_fraction:.2f}. "
                f"The chosen option has a fraction of {chosen_fraction}."
            )
        if not is_power_correct:
            error_messages.append(
                f"Constraint violated: Wavelength dependence. "
                f"The power should be proportional to λ^({expected_power_of_lambda}). "
                f"The chosen option has a dependence of λ^({chosen_power})."
            )
        return "\n".join(error_messages)

# Run the check
result = check_dipole_radiation_answer()
print(result)