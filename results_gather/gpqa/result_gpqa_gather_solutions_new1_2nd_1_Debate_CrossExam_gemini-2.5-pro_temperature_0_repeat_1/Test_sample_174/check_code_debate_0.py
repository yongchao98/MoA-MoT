import math
import re

def check_correctness_of_answer():
    """
    Checks the correctness of the final answer based on physics principles.
    
    An oscillating spheroid is modeled as an electric dipole.
    1.  The radiated power per unit solid angle is proportional to sin^2(theta).
        The maximum is at theta=90 (sin^2(90)=1).
        At theta=30, the fraction is sin^2(30) = (1/2)^2 = 1/4.
    2.  The radiated power is proportional to omega^4, and since omega is proportional
        to 1/lambda, the power is proportional to lambda^-4.
    
    The correct answer should therefore be (1/4, lambda^-4).
    """
    
    # --- Step 1: Calculate the correct values from physics principles ---
    
    # Calculate the correct fraction
    theta_degrees = 30
    theta_radians = math.radians(theta_degrees)
    # The fraction of maximum power is sin^2(theta) / sin^2(90)
    correct_fraction = math.sin(theta_radians)**2
    
    # Determine the correct power of the wavelength
    # Power is proportional to omega^4, and omega is proportional to 1/lambda.
    # So, Power is proportional to (1/lambda)^4 = lambda^-4.
    correct_wavelength_power = -4

    # --- Step 2: Parse the provided options and the final answer ---
    
    # The options as listed in the question
    # Note: λ^(-4) is the same as λ⁻⁴
    options = {
        "A": {"fraction": 1/2, "power": -4},
        "B": {"fraction": 1/4, "power": -4},
        "C": {"fraction": 3/4, "power": -6},
        "D": {"fraction": 1/4, "power": -3}
    }
    
    # The final answer provided in the prompt
    final_answer_str = "<<<B>>>"
    
    # Extract the letter from the final answer string
    match = re.search(r'<<<([A-D])>>>', final_answer_str)
    if not match:
        return "Error: Could not parse the final answer format."
        
    chosen_option_letter = match.group(1)
    
    if chosen_option_letter not in options:
        return f"Error: The chosen option '{chosen_option_letter}' is not a valid option."

    chosen_option_values = options[chosen_option_letter]
    
    # --- Step 3: Compare the calculated correct answer with the chosen option ---
    
    # Check if the fraction matches
    if not math.isclose(correct_fraction, chosen_option_values["fraction"]):
        return (f"Incorrect: The fraction of power is wrong. "
                f"Based on physics (sin^2(30°)), the expected fraction is {correct_fraction:.2f}. "
                f"The chosen option '{chosen_option_letter}' provides a fraction of {chosen_option_values['fraction']}.")

    # Check if the wavelength power dependence matches
    if correct_wavelength_power != chosen_option_values["power"]:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"Based on physics (Power ∝ λ^n), the expected power 'n' is {correct_wavelength_power}. "
                f"The chosen option '{chosen_option_letter}' provides a power of {chosen_option_values['power']}.")

    # If both checks pass, the answer is correct
    return "Correct"

# Run the check and print the result
result = check_correctness_of_answer()
print(result)