import numpy as np

def check_correctness_of_llm_answer():
    """
    This function checks the correctness of the LLM's answer to the physics problem.
    The problem asks for the fraction of maximum radiated power at 30 degrees and the
    wavelength dependence for an oscillating spheroid.
    The LLM's answer is 'A', which corresponds to a fraction of 1/4 and a wavelength
    dependence of lambda^(-4).
    """
    
    # The LLM's chosen answer is A.
    llm_answer_option = "A"
    
    # Define the options from the question.
    options = {
        "A": {"fraction": 1/4, "lambda_power": -4},
        "B": {"fraction": 3/4, "lambda_power": -6},
        "C": {"fraction": 1/4, "lambda_power": -3},
        "D": {"fraction": 1/2, "lambda_power": -4},
    }
    
    # Get the values corresponding to the LLM's answer.
    llm_answer_values = options[llm_answer_option]
    llm_fraction = llm_answer_values["fraction"]
    llm_lambda_power = llm_answer_values["lambda_power"]

    # --- THEORETICAL CALCULATION BASED ON PHYSICS ---
    # The most general oscillation of a charge distribution will radiate primarily
    # as an electric dipole, which is the lowest-order term.

    # 1. Calculate the theoretical fraction of maximum power.
    # For an electric dipole oscillating along the z-axis, the radiated power
    # per unit solid angle is dP/dΩ ∝ sin²(θ).
    # The maximum power occurs at θ = 90° (where sin²(90°) = 1).
    # The fraction is the ratio of the power at 30° to the maximum power.
    theta_rad = np.deg2rad(30)
    theta_max_power_rad = np.deg2rad(90)
    
    theoretical_fraction = (np.sin(theta_rad)**2) / (np.sin(theta_max_power_rad)**2)

    # 2. Determine the theoretical wavelength dependence.
    # For electric dipole radiation, the total radiated power P is proportional
    # to ω⁴. Since ω = 2πc/λ, the power is proportional to (1/λ)⁴ = λ⁻⁴.
    theoretical_lambda_power = -4

    # --- VERIFICATION ---
    # Check if the LLM's answer matches the theoretical values.
    
    # Check the fraction part.
    if not np.isclose(llm_fraction, theoretical_fraction):
        return (f"Incorrect: The fraction of maximum power is wrong. "
                f"The LLM's answer states the fraction is {llm_fraction}. "
                f"However, for the dominant electric dipole radiation model, the power is proportional to sin²(θ). "
                f"The theoretical fraction at θ=30° is sin²(30°)/sin²(90°) = {theoretical_fraction:.4f}. "
                f"The answer's fraction does not satisfy this constraint.")

    # Check the wavelength dependence part.
    if llm_lambda_power != theoretical_lambda_power:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"The LLM's answer states the power is proportional to λ^({llm_lambda_power}). "
                f"However, for electric dipole radiation, the power is proportional to λ^({theoretical_lambda_power}). "
                f"The answer's wavelength dependence does not satisfy this constraint.")

    # If both parts are correct, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_correctness_of_llm_answer()
print(result)