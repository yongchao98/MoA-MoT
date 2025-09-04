import math

def check_answer_correctness():
    """
    Checks the correctness of the answer based on the physics of electric dipole radiation.
    
    The function verifies two key aspects:
    1. The fraction of maximum power radiated at a specific angle (30 degrees).
    2. The dependence of radiated power on the wavelength (lambda).
    """

    # --- Step 1: Define the options from the question ---
    options = {
        "A": {"fraction": 1/4, "lambda_exp": -3},
        "B": {"fraction": 1/2, "lambda_exp": -4},
        "C": {"fraction": 1/4, "lambda_exp": -4},
        "D": {"fraction": 3/4, "lambda_exp": -6}
    }
    
    # The final answer provided by the LLM to be checked.
    llm_answer = "C"

    # --- Step 2: Calculate the fraction based on the angular distribution ---
    # For electric dipole radiation, power is proportional to sin^2(theta).
    # Maximum power is at theta = 90 degrees.
    theta_max_rad = math.radians(90)
    # Power is measured at theta = 30 degrees.
    theta_test_rad = math.radians(30)

    # The proportional power values are sin^2(theta).
    power_max_proportional = math.sin(theta_max_rad)**2
    power_test_proportional = math.sin(theta_test_rad)**2

    # The fraction is the ratio of the power at 30 degrees to the maximum power.
    calculated_fraction = power_test_proportional / power_max_proportional

    # --- Step 3: Determine the wavelength dependence exponent ---
    # For electric dipole radiation, power is proportional to omega^4.
    # Since omega is proportional to 1/lambda, power is proportional to (1/lambda)^4 = lambda^-4.
    calculated_lambda_exponent = -4

    # --- Step 4: Compare the calculated results with the chosen option ---
    chosen_option_values = options.get(llm_answer)
    
    if not chosen_option_values:
        return f"Invalid answer key '{llm_answer}' provided. It is not one of the options A, B, C, D."

    expected_fraction = chosen_option_values["fraction"]
    expected_lambda_exp = chosen_option_values["lambda_exp"]

    # Check if both parts of the answer are correct
    is_fraction_correct = math.isclose(calculated_fraction, expected_fraction)
    is_lambda_exp_correct = (calculated_lambda_exponent == expected_lambda_exp)

    if is_fraction_correct and is_lambda_exp_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_fraction_correct:
            error_messages.append(
                f"The fraction of maximum power is incorrect. "
                f"Based on the sin^2(theta) dependence, the fraction at 30 degrees should be {calculated_fraction:.2f}, "
                f"but option {llm_answer} states it is {expected_fraction:.2f}."
            )
        if not is_lambda_exp_correct:
            error_messages.append(
                f"The wavelength dependence is incorrect. "
                f"For electric dipole radiation, power is proportional to lambda^{calculated_lambda_exponent}, "
                f"but option {llm_answer} states it is lambda^{expected_lambda_exp}."
            )
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_answer_correctness()
print(result)