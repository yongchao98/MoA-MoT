import math

def check_exoplanet_period_ratio():
    """
    Checks the correctness of the answer for the exoplanet period ratio problem.

    The function performs the following steps:
    1. Defines the given parameters: wavelength shifts and multiple-choice options.
    2. Calculates the theoretical ratio of the orbital periods (P₂ / P₁) based on the physics of the radial velocity method.
       - The radial velocity semi-amplitude K is proportional to the wavelength shift Δλ.
       - K is also proportional to P^(-1/3), where P is the orbital period.
       - This leads to the relationship: P₂ / P₁ = (Δλ₁ / Δλ₂)³.
    3. Compares the calculated result with the value of the provided answer ('D').
    4. Verifies that the provided answer is the closest match among all options.
    """
    # Given parameters from the question
    delta_lambda_1 = 5.0  # miliangstrom
    delta_lambda_2 = 7.0  # miliangstrom
    
    # Multiple-choice options provided in the question
    options = {
        'A': 1.96,
        'B': 0.85,
        'C': 1.40,
        'D': 0.36
    }
    
    # The final answer to be checked
    provided_answer_letter = 'D'

    # Step 1: Calculate the theoretical ratio of the orbital periods
    # K₂ / K₁ = Δλ₂ / Δλ₁
    # K₂ / K₁ = (P₂ / P₁)^(-1/3)
    # (Δλ₂ / Δλ₁)^(-3) = P₂ / P₁
    # P₂ / P₁ = (Δλ₁ / Δλ₂)^3
    try:
        calculated_period_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. delta_lambda_2 cannot be zero."
    
    # Step 2: Check if the provided answer is a valid option
    if provided_answer_letter not in options:
        return f"The provided answer '{provided_answer_letter}' is not one of the valid options {list(options.keys())}."

    provided_answer_value = options[provided_answer_letter]

    # Step 3: Check if the calculated value matches the provided answer's value
    # We use a tolerance to account for potential rounding in the options
    tolerance = 0.01 
    if not math.isclose(calculated_period_ratio, provided_answer_value, abs_tol=tolerance):
        return (f"The calculated period ratio is {calculated_period_ratio:.4f}. "
                f"The value for the provided answer '{provided_answer_letter}' is {provided_answer_value}, which is not a close match.")

    # Step 4: Verify that the provided answer is the BEST match among all options
    best_fit_option = min(options, key=lambda k: abs(options[k] - calculated_period_ratio))

    if best_fit_option != provided_answer_letter:
        return (f"The calculated period ratio is {calculated_period_ratio:.4f}. "
                f"The best matching option is '{best_fit_option}' (value: {options[best_fit_option]}), "
                f"but the provided answer was '{provided_answer_letter}'.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_exoplanet_period_ratio()
print(result)