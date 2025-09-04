import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.
    
    The problem asks for the angular distance between the first two minima for
    diffraction from a circular aperture of radius 'a'. The formula for this is:
    Δθ = ((z₂ - z₁) / 2π) * (λ/a)
    where z₁ and z₂ are the first two non-trivial zeros of the J₁(x) Bessel function.
    """
    
    # --- Step 1: Define the problem parameters from the question ---
    
    # Options provided in the question text
    # A) 0.610 λ / a
    # B) 0.506 λ / a
    # C) 1.220 λ /a
    # D) 0.500 λ / a
    options = {
        'A': 0.610,
        'B': 0.506,
        'C': 1.220,
        'D': 0.500
    }
    
    # The final answer given by the LLM is <<<B>>>
    proposed_answer_key = 'B'
    
    # --- Step 2: Perform the correct physical calculation ---
    
    try:
        # jn_zeros(n, k) finds the first k positive zeros of the Bessel function Jn(x).
        # We need the first 2 zeros of J1(x).
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]  # First zero
        z2 = zeros[1]  # Second zero
    except (ImportError, ModuleNotFoundError):
        # In case scipy is not available, use well-known values.
        print("SciPy not found, using hardcoded values for Bessel zeros.")
        z1 = 3.83170597
        z2 = 7.01558667

    # Calculate the theoretical coefficient for the angular distance Δθ
    theoretical_coefficient = (z2 - z1) / (2 * np.pi)
    
    # --- Step 3: Verify the proposed answer ---
    
    # Check if the proposed answer key is valid
    if proposed_answer_key not in options:
        return f"Invalid answer key '{proposed_answer_key}'. Valid keys are {list(options.keys())}."

    proposed_coefficient = options[proposed_answer_key]

    # Find the option that is numerically closest to the theoretical value
    best_match_key = min(options, key=lambda k: abs(options[k] - theoretical_coefficient))

    # Check for common mistakes to provide a more detailed reason for failure
    # Mistake 1: Reporting the position of the first minimum (θ₁)
    theta1_coeff = z1 / (2 * np.pi)
    
    # Mistake 2: Reporting the angular diameter of the central maximum (2θ₁)
    diameter_coeff = 2 * z1 / (2 * np.pi)

    if proposed_answer_key == best_match_key:
        # The reasoning in the provided text is also sound. It correctly identifies
        # the geometry, the physics (Fraunhofer diffraction, Bessel zeros), and the
        # quantity to be calculated (θ₂ - θ₁).
        return "Correct"
    else:
        reason = (
            f"The proposed answer is '{proposed_answer_key}' ({proposed_coefficient}), but the correct answer is '{best_match_key}' ({options[best_match_key]}).\n"
            f"The theoretical coefficient for the angular distance is (z₂ - z₁) / 2π ≈ {theoretical_coefficient:.5f}.\n"
            f"This value is closest to {options[best_match_key]}."
        )
        if abs(proposed_coefficient - theta1_coeff) < 0.001:
            reason += f"\nThe proposed answer '{proposed_answer_key}' ({proposed_coefficient}) corresponds to the position of the first minimum (θ₁ ≈ {theta1_coeff:.3f} λ/a), not the angular distance between the first two minima."
        elif abs(proposed_coefficient - diameter_coeff) < 0.001:
            reason += f"\nThe proposed answer '{proposed_answer_key}' ({proposed_coefficient}) corresponds to the angular diameter of the central maximum (2θ₁ ≈ {diameter_coeff:.3f} λ/a), not the angular distance between the first two minima."
        return reason

# Execute the check and print the result
result = check_diffraction_answer()
print(result)