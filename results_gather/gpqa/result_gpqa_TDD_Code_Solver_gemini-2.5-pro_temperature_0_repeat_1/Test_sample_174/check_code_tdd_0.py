import math

def check_spheroid_radiation_answer():
    """
    Checks the correctness of the provided answer for the oscillating spheroid problem.
    """
    # --- Problem Parameters ---
    theta_degrees = 30.0
    llm_answer_key = "D"

    # --- Physics Principles ---
    # 1. Calculate the expected fraction of maximum power.
    # The power is proportional to sin^2(theta). Maximum power is at theta=90 deg.
    # Fraction = sin^2(theta_given) / sin^2(90) = sin^2(theta_given)
    theta_radians = math.radians(theta_degrees)
    expected_fraction = math.sin(theta_radians)**2

    # 2. Determine the expected wavelength dependence.
    # Power is proportional to omega^4, and omega is proportional to 1/lambda.
    # So, Power is proportional to (1/lambda)^4 = lambda^-4.
    expected_lambda_exponent = -4

    # --- Options from the question ---
    options = {
        "A": {"fraction": 1/4, "lambda_exponent": -3},
        "B": {"fraction": 3/4, "lambda_exponent": -6},
        "C": {"fraction": 1/2, "lambda_exponent": -4},
        "D": {"fraction": 1/4, "lambda_exponent": -4}
    }

    # --- Verification ---
    if llm_answer_key not in options:
        return f"Invalid answer key '{llm_answer_key}'. The key must be one of {list(options.keys())}."

    chosen_option = options[llm_answer_key]

    # Check if the fraction in the chosen option is correct
    if not math.isclose(chosen_option["fraction"], expected_fraction, rel_tol=1e-9):
        return (f"Incorrect: The fraction of maximum power is wrong. "
                f"For θ = 30°, the fraction should be sin²(30°) = (0.5)² = 0.25 (or 1/4). "
                f"The chosen option '{llm_answer_key}' has a fraction of {chosen_option['fraction']}.")

    # Check if the lambda exponent in the chosen option is correct
    if chosen_option["lambda_exponent"] != expected_lambda_exponent:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"The radiated power is proportional to λ⁻⁴, so the exponent should be -4. "
                f"The chosen option '{llm_answer_key}' has an exponent of {chosen_option['lambda_exponent']}.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_spheroid_radiation_answer()
print(result)