import numpy as np
try:
    from scipy.special import jn_zeros
    scipy_available = True
except ImportError:
    scipy_available = False

def check_diffraction_answer():
    """
    Checks the correctness of the provided answer for the physics problem.

    The problem asks for the coefficient C in the expression for the angular
    distance between the first two diffraction minima: Δθ = C * (λ / a).
    """
    # The final answer provided is 'D', which corresponds to the coefficient 0.506.
    answer_coefficient = 0.506

    # --- Theoretical Calculation ---
    # The coefficient C is given by (z₂ - z₁) / (2π), where z₁ and z₂
    # are the first two non-trivial zeros of the J₁ Bessel function.

    if scipy_available:
        # Use scipy for high-precision values of the Bessel zeros.
        # jn_zeros(1, 2) computes the first 2 positive zeros of J₁(x).
        zeros = jn_zeros(1, 2)
        z1, z2 = zeros[0], zeros[1]
    else:
        # Fallback to well-known approximate values if scipy is not installed.
        print("Warning: SciPy not found. Using approximate values for Bessel zeros.")
        z1 = 3.83170597
        z2 = 7.01558667

    # Calculate the theoretical coefficient.
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # --- Verification ---
    # The options are given to 3 decimal places. We check if the answer's
    # coefficient is the closest option to the calculated theoretical value.
    options = {
        'A': 0.500,
        'B': 1.220,
        'C': 0.610,
        'D': 0.506
    }

    # Find the option with the minimum difference from the calculated value.
    best_option = min(options, key=lambda k: abs(options[k] - calculated_coefficient))

    # The provided answer is 'D'. Check if our calculation agrees.
    if best_option == 'D':
        return "Correct"
    else:
        reason = (
            f"The provided answer 'D' (coefficient {answer_coefficient}) is incorrect.\n"
            f"The theoretical calculation is as follows:\n"
            f"1. The aperture is a circle of radius 'a'.\n"
            f"2. The angular distance between the first two minima is Δθ = [(z₂ - z₁) / (2π)] * (λ / a).\n"
            f"3. The first two zeros of the J₁ Bessel function are z₁ ≈ {z1:.6f} and z₂ ≈ {z2:.6f}.\n"
            f"4. The calculated coefficient is C = ({z2:.6f} - {z1:.6f}) / (2π) ≈ {calculated_coefficient:.6f}.\n"
            f"5. The calculated coefficient {calculated_coefficient:.6f} is closest to option '{best_option}' (value: {options[best_option]}), not 'D'."
        )
        return reason

# Run the check and print the result.
print(check_diffraction_answer())