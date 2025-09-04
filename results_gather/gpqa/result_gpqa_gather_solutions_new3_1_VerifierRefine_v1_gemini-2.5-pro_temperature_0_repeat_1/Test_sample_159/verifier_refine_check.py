import math

def check_answer():
    """
    Checks the correctness of the answer for the diffraction problem.
    """
    # The final answer from the LLM is 'A', which corresponds to the coefficient 0.506.
    # Options are: A) 0.506, B) 0.610, C) 0.500, D) 1.220
    # The provided answer is 'A'.
    expected_coefficient = 0.506
    
    # We use accurate values for the first two non-trivial zeros of the 
    # first-order Bessel function J1(x).
    # z1 corresponds to the first minimum, z2 to the second.
    try:
        # Use scipy for the most accurate values if available
        from scipy.special import jn_zeros
        # jn_zeros(n, k) gives the first k positive zeros of Jn(x)
        zeros = jn_zeros(1, 2)
        z1, z2 = zeros[0], zeros[1]
    except ImportError:
        # Fallback to high-precision hardcoded values if scipy is not installed
        z1 = 3.8317059702075123
        z2 = 7.015586669815613

    # The angular distance between the first two minima is given by:
    # Δθ = θ₂ - θ₁ = ((z₂ - z₁) / (2π)) * (λ/a)
    # We need to calculate the coefficient C = (z₂ - z₁) / (2π)
    calculated_coefficient = (z2 - z1) / (2 * math.pi)

    # Check if the calculated coefficient matches the expected one from option A.
    # The options are given to 3 decimal places, so a tolerance of 0.001 is appropriate.
    tolerance = 0.001
    
    if abs(calculated_coefficient - expected_coefficient) < tolerance:
        return "Correct"
    else:
        # Check other options to provide a more detailed reason for failure.
        first_min_coeff = z1 / (2 * math.pi) # This is approx 0.610 (Option B)
        
        if abs(first_min_coeff - expected_coefficient) < tolerance:
             return (f"Incorrect. The provided answer's coefficient {expected_coefficient} corresponds to the "
                    f"angular position of the *first* minimum (θ₁ ≈ {first_min_coeff:.4f} λ/a), "
                    f"not the angular distance between the first two minima (Δθ = θ₂ - θ₁).")

        return (f"Incorrect. The calculated coefficient for the angular distance is "
                f"C = (z₂ - z₁) / (2π) ≈ {calculated_coefficient:.4f}. The coefficient from the "
                f"selected answer is {expected_coefficient}. These values do not match within the "
                f"required tolerance.")

# Run the check
result = check_answer()
print(result)