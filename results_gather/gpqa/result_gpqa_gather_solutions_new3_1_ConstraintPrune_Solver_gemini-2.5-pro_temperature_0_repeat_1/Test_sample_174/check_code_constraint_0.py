import math

def check_answer():
    """
    Checks the correctness of the answer based on the physics of electric dipole radiation.
    """
    # The final answer to be checked.
    # The provided solution identifies 'D' as the correct option.
    final_answer_key = "D"

    # Define the options from the question
    options = {
        "A": {"fraction": 3/4, "lambda_power": -6},
        "B": {"fraction": 1/4, "lambda_power": -3},
        "C": {"fraction": 1/2, "lambda_power": -4},
        "D": {"fraction": 1/4, "lambda_power": -4},
    }

    # --- Constraint 1: Calculate the expected fraction of power at 30 degrees ---
    # Power is proportional to sin^2(theta).
    # Maximum power is at theta = 90 degrees, where sin^2(90) = 1.
    # The fraction at 30 degrees is sin^2(30) / sin^2(90).
    theta_rad = math.radians(30)
    expected_fraction = math.sin(theta_rad)**2

    # --- Constraint 2: Define the expected wavelength dependence ---
    # Power is proportional to omega^4, and omega is proportional to 1/lambda.
    # Therefore, power is proportional to (1/lambda)^4 = lambda^-4.
    expected_lambda_power = -4

    # Get the data for the proposed answer
    proposed_answer_data = options.get(final_answer_key)

    if not proposed_answer_data:
        return f"Error: The answer key '{final_answer_key}' is not a valid option."

    # Check if the proposed answer satisfies the constraints
    fraction_matches = math.isclose(proposed_answer_data["fraction"], expected_fraction)
    lambda_power_matches = (proposed_answer_data["lambda_power"] == expected_lambda_power)

    # Formulate the result
    if fraction_matches and lambda_power_matches:
        return "Correct"
    else:
        error_reasons = []
        if not fraction_matches:
            error_reasons.append(
                f"the fraction of power at 30 degrees is incorrect. "
                f"Expected {expected_fraction:.2f} (1/4), but option {final_answer_key} gives {proposed_answer_data['fraction']}."
            )
        if not lambda_power_matches:
            error_reasons.append(
                f"the wavelength dependence is incorrect. "
                f"Expected a form of lambda^{expected_lambda_power}, but option {final_answer_key} gives lambda^{proposed_answer_data['lambda_power']}."
            )
        return "Incorrect because " + " and ".join(error_reasons)

# Run the check and print the result
result = check_answer()
print(result)