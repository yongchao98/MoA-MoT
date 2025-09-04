import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.

    The problem asks for the angular distance between the first two minima for diffraction
    through a circular aperture of radius 'a'. The formula for this distance is:
    Δθ = ((z₂ - z₁) / (2π)) * (λ / a)
    where z₁ and z₂ are the first two non-trivial zeros of the first-order Bessel function J₁(x).

    This function calculates the theoretical coefficient C = (z₂ - z₁) / (2π) and compares it
    to the options provided in the question.
    """
    # Step 1: Define the options from the question prompt.
    # Note: The final answer text has a different ordering of options than some of the
    # candidate answers. We must use the ordering from the final answer's analysis.
    # A) 0.610, B) 0.500, C) 1.220, D) 0.506
    options = {
        "A": 0.610,
        "B": 0.500,
        "C": 1.220,
        "D": 0.506
    }
    
    # The final answer selected is <<<D>>>.
    selected_option_letter = "D"
    answer_coeff = options.get(selected_option_letter)

    # Step 2: Calculate the theoretical coefficient using SciPy.
    try:
        # jn_zeros(n, k) returns the first k positive zeros of the Bessel function Jn(x).
        # We need the first 2 zeros of J1(x).
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]
        z2 = zeros[1]
    except ImportError:
        return "Scipy library not found. Cannot perform the check."
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    theoretical_coeff = (z2 - z1) / (2 * np.pi)

    # Step 3: Verify the logic and correctness.
    # The analysis in the provided answer is correct.
    # 1. It correctly identifies that an N-sided polygon with constant apothem 'a'
    #    becomes a circle of radius 'a'.
    # 2. It correctly identifies the formula for the minima based on the zeros of the
    #    Bessel function J1(x).
    # 3. It correctly identifies that the question asks for the difference between the
    #    first two minima (θ₂ - θ₁).
    # 4. The calculation yields ~0.5068.

    # Step 4: Check if the selected option is the best fit.
    # Since the options are rounded, we find which option is closest to the theoretical value.
    differences = {key: abs(value - theoretical_coeff) for key, value in options.items()}
    best_option_letter = min(differences, key=differences.get)

    if best_option_letter == selected_option_letter:
        # The selected option is indeed the closest to the theoretical value.
        # The difference between the theoretical value (~0.50673) and the chosen answer (0.506)
        # is ~0.00073, which is a reasonable rounding difference.
        return "Correct"
    else:
        return (f"Incorrect. The theoretical coefficient is approximately {theoretical_coeff:.5f}. "
                f"The closest option is '{best_option_letter}' with a value of {options[best_option_letter]}. "
                f"The provided answer selected '{selected_option_letter}' with a value of {answer_coeff}.")

# Execute the check and print the result.
result = check_diffraction_answer()
print(result)