import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the angular distance between the first two minima
    in the diffraction pattern of a circular aperture.
    """
    try:
        # The problem asks for the coefficient C in the expression:
        # angular_distance = C * (lambda / a)
        # where C = (z2 - z1) / (2 * pi)
        # z1 and z2 are the first two non-trivial zeros of the Bessel function J1(x).

        # Get the first two non-trivial zeros of J1(x) using SciPy
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]  # approx 3.8317
        z2 = zeros[1]  # approx 7.0156

        # Calculate the correct coefficient for the angular distance
        correct_coefficient = (z2 - z1) / (2 * np.pi)

        # The options provided in the question
        options = {
            "A": 0.500,
            "B": 1.220,
            "C": 0.506,
            "D": 0.610
        }
        
        # The final answer provided by the LLM is C
        llm_answer_choice = "C"
        llm_answer_value = options[llm_answer_choice]

        # Find the option that is numerically closest to the calculated correct value
        best_option_choice = min(options, key=lambda k: abs(options[k] - correct_coefficient))

        # Check if the LLM's answer matches the best option
        if llm_answer_choice == best_option_choice:
            # The LLM chose the correct option.
            # Let's verify the precision. The difference should be small.
            if abs(llm_answer_value - correct_coefficient) < 0.001:
                return "Correct"
            else:
                # This case is unlikely given the options
                return f"The chosen answer {llm_answer_choice} ({llm_answer_value}) is the closest option, but its value is not precise. The calculated value is {correct_coefficient:.5f}."
        else:
            # The LLM chose the wrong option. Provide a reason.
            # Check for common mistakes.
            
            # Mistake 1: Reporting the position of the first minimum, theta_1
            first_min_coeff = z1 / (2 * np.pi) # This should be ~0.610
            if abs(llm_answer_value - first_min_coeff) < 0.001:
                return f"Incorrect. The answer {llm_answer_choice} ({llm_answer_value}) corresponds to the angular position of the *first* minimum (approx {first_min_coeff:.3f} * lambda/a), not the angular distance between the first two minima."

            # Mistake 2: Using the coefficient for the diameter D (1.220)
            # theta_1 = 1.220 * lambda / D = 1.220 * lambda / (2a) = 0.610 * lambda / a
            # So, a coefficient of 1.220 for lambda/a is incorrect.
            if llm_answer_choice == "B":
                 return f"Incorrect. The value 1.220 is the coefficient for the first minimum when using the diameter D (i.e., theta_1 = 1.220 * lambda/D). The question asks for a coefficient for lambda/a, and for the distance between the first two minima."

            return f"Incorrect. The chosen answer {llm_answer_choice} ({llm_answer_value}) is wrong. The correct coefficient is ({z2:.4f} - {z1:.4f}) / (2 * pi) = {correct_coefficient:.5f}, which is best matched by option {best_option_choice} ({options[best_option_choice]})."

    except ImportError:
        return "Verification failed: The 'scipy' library is required. Please install it using 'pip install scipy numpy'."
    except Exception as e:
        return f"An error occurred during verification: {e}"

# Run the check
result = check_diffraction_answer()
print(result)