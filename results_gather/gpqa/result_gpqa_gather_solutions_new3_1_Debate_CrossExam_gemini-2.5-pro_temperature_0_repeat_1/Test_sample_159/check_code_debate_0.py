import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.

    The problem asks for the angular distance between the first two minima for a
    circular aperture of radius 'a'. The formula for this distance is:
    Δθ = (z₂ - z₁) / (2π) * (λ / a)
    where z₁ and z₂ are the first two non-trivial zeros of the Bessel function J₁(x).

    This function calculates the coefficient C = (z₂ - z₁) / (2π) and compares it
    to the coefficient from the provided answer.
    """
    
    # The final answer provided by the LLM is 'B', which corresponds to a coefficient of 0.506.
    llm_answer_label = 'B'
    options = {
        'A': 1.220,
        'B': 0.506,
        'C': 0.610,
        'D': 0.500
    }
    
    if llm_answer_label not in options:
        return f"Invalid answer label '{llm_answer_label}'. The label must be one of {list(options.keys())}."
        
    llm_answer_coefficient = options[llm_answer_label]

    # Step 1: Find the first two non-trivial zeros of the Bessel function J₁(x).
    # jn_zeros(1, 2) finds the first 2 zeros of the Bessel function of order 1.
    try:
        zeros = jn_zeros(1, 2)
        z1, z2 = zeros[0], zeros[1]
    except ImportError:
        # Fallback to high-precision known values if scipy is not available.
        z1 = 3.8317059702075123
        z2 = 7.015586669815618
    
    # Step 2: Calculate the theoretical coefficient for the angular distance.
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # Step 3: Check if the LLM's answer is the closest one to the calculated value.
    # Find which option is numerically closest to our calculation.
    closest_label = min(options, key=lambda k: abs(options[k] - calculated_coefficient))

    if closest_label == llm_answer_label:
        # The LLM correctly identified the closest option.
        return "Correct"
    else:
        # The LLM chose the wrong option.
        # Let's provide a detailed reason.
        
        # Check for common mistakes to provide a better explanation.
        # Mistake 1: Position of the first minimum, θ₁ = (z₁ / 2π) * (λ/a)
        first_min_coeff = z1 / (2 * np.pi)
        if np.isclose(llm_answer_coefficient, first_min_coeff, atol=0.001):
            return (f"Incorrect. The provided answer ({llm_answer_coefficient}) corresponds to the angular position of the "
                    f"first minimum (approx. {first_min_coeff:.3f} λ/a), not the angular distance between the first two minima. "
                    f"The correct answer is approx. {calculated_coefficient:.3f} λ/a.")

        # Mistake 2: Angular diameter of the central maximum, 2θ₁ = (z₁ / π) * (λ/a)
        central_max_diam_coeff = z1 / np.pi
        if np.isclose(llm_answer_coefficient, central_max_diam_coeff, atol=0.001):
            return (f"Incorrect. The provided answer ({llm_answer_coefficient}) corresponds to the angular diameter of the "
                    f"central maximum (approx. {central_max_diam_coeff:.3f} λ/a), not the angular distance between the first two minima. "
                    f"The correct answer is approx. {calculated_coefficient:.3f} λ/a.")

        return (f"Incorrect. The calculated coefficient for the angular distance is (z₂ - z₁) / (2π) ≈ {calculated_coefficient:.4f}. "
                f"Among the given options, this value is closest to option '{closest_label}' ({options[closest_label]}), "
                f"not the provided answer '{llm_answer_label}' ({llm_answer_coefficient}).")

# Run the check
result = check_diffraction_answer()
print(result)