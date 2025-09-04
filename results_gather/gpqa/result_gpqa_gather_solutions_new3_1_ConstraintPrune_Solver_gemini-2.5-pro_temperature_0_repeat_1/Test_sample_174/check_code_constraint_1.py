import math

def check_correctness():
    """
    Checks the correctness of the answer based on the physics of electric dipole radiation.
    """
    # The final answer provided by the LLM to be checked.
    llm_final_answer = "D"

    # Define the options from the question.
    options = {
        "A": {"fraction": 3/4, "lambda_exponent": -6},
        "B": {"fraction": 1/4, "lambda_exponent": -3},
        "C": {"fraction": 1/2, "lambda_exponent": -4},
        "D": {"fraction": 1/4, "lambda_exponent": -4}
    }

    # --- Constraint 1: Calculate the correct fraction for the angular dependence ---
    # The radiated power is proportional to sin^2(theta).
    # The maximum is at theta=90 degrees (sin^2(90)=1).
    # We need the fraction at theta=30 degrees.
    angle_deg = 30
    angle_rad = math.radians(angle_deg)
    
    # The fraction of maximum power is sin^2(30) / sin^2(90)
    correct_fraction = math.sin(angle_rad)**2

    # --- Constraint 2: Determine the correct wavelength dependence ---
    # Radiated power is proportional to omega^4.
    # Since omega is proportional to 1/lambda, power is proportional to lambda^-4.
    correct_lambda_exponent = -4

    # --- Verification ---
    # Get the values from the LLM's chosen option.
    chosen_option_data = options.get(llm_final_answer)

    if not chosen_option_data:
        return f"Incorrect. The provided answer '{llm_final_answer}' is not a valid option key."

    # Check if the fraction matches the calculated correct fraction.
    if not math.isclose(chosen_option_data["fraction"], correct_fraction):
        return (f"Incorrect. The fraction of maximum power is wrong. "
                f"The correct fraction for θ=30° is sin²(30°) = {correct_fraction:.2f}. "
                f"The answer's option '{llm_final_answer}' provides a fraction of {chosen_option_data['fraction']}.")

    # Check if the lambda exponent matches the correct exponent.
    if chosen_option_data["lambda_exponent"] != correct_lambda_exponent:
        return (f"Incorrect. The wavelength dependence is wrong. "
                f"The radiated power for an electric dipole is proportional to λ^({correct_lambda_exponent}). "
                f"The answer's option '{llm_final_answer}' provides a dependence of λ^({chosen_option_data['lambda_exponent']}).")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)