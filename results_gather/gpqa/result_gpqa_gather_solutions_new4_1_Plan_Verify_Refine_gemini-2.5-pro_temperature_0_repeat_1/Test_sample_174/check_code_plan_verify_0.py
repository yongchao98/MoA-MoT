import math

def check_correctness():
    """
    Checks the correctness of the answer based on the principles of electric dipole radiation.
    """
    # The final answer provided by the LLM analysis.
    final_answer = "A"

    # --- Step 1: Define the options from the question ---
    options = {
        "A": {"fraction": 1/4, "lambda_exp": -4},
        "B": {"fraction": 1/2, "lambda_exp": -4},
        "C": {"fraction": 1/4, "lambda_exp": -3},
        "D": {"fraction": 3/4, "lambda_exp": -6}
    }

    # --- Step 2: Calculate the expected values based on physics principles ---
    # The standard model for an oscillating charge distribution is an electric dipole.

    # Part a: Calculate the fraction of maximum power at theta = 30 degrees.
    # Power is proportional to sin^2(theta).
    # Maximum power is at theta = 90 degrees, where sin^2(90) = 1.
    # The fraction is sin^2(30) / sin^2(90) = sin^2(30).
    theta_deg = 30
    theta_rad = math.radians(theta_deg)
    expected_fraction = math.sin(theta_rad)**2

    # Part b: Determine the wavelength dependence.
    # Radiated power is proportional to omega^4.
    # Since omega is proportional to 1/lambda, power is proportional to lambda^-4.
    expected_lambda_exp = -4

    # --- Step 3: Compare the LLM's answer with the calculated values ---
    chosen_option_data = options.get(final_answer)

    if not chosen_option_data:
        return f"Invalid answer format. The answer '{final_answer}' is not one of the options A, B, C, D."

    # Check if the fraction in the chosen option is correct.
    if not math.isclose(chosen_option_data["fraction"], expected_fraction):
        return (f"Incorrect. The fraction of maximum power is wrong. "
                f"For an electric dipole, power is proportional to sin^2(theta). "
                f"The maximum is at theta=90 deg (sin^2(90)=1). "
                f"At theta=30 deg, the fraction of maximum power is sin^2(30) = (1/2)^2 = 0.25. "
                f"The answer '{final_answer}' states the fraction is {chosen_option_data['fraction']}, but it should be {expected_fraction:.2f}.")

    # Check if the lambda dependence in the chosen option is correct.
    if chosen_option_data["lambda_exp"] != expected_lambda_exp:
        return (f"Incorrect. The wavelength dependence is wrong. "
                f"For an electric dipole, radiated power is proportional to omega^4, "
                f"which means it is proportional to lambda^({expected_lambda_exp}). "
                f"The answer '{final_answer}' states the dependence is lambda^({chosen_option_data['lambda_exp']}).")

    # If both parts of the answer are correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)