import math

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the angular distance between the first two minima
    in a circular aperture diffraction problem.
    """
    
    # The final answer provided is 'C', which corresponds to a coefficient of 0.506.
    answer_coefficient = 0.506
    
    # --- Theoretical Calculation ---
    # The angular distance is given by Δθ = (z₂ - z₁) / (2π) * (λ / a).
    # We need the first two non-zero roots (zeros) of the Bessel function J₁(x).
    
    # Use high-precision values for the zeros.
    # For a more robust check, one could use scipy.special.jn_zeros(1, 2),
    # but for self-containment, we use well-known values.
    z1 = 3.831706  # First non-zero root of J₁(x)
    z2 = 7.015587  # Second non-zero root of J₁(x)
    
    # Calculate the theoretical coefficient C = (z₂ - z₁) / (2π)
    theoretical_coefficient = (z2 - z1) / (2 * math.pi)
    
    # --- Verification ---
    # Compare the theoretical coefficient with the one from the chosen answer 'C'.
    # We use math.isclose() to handle potential floating-point inaccuracies and rounding in the option.
    # A relative tolerance of 1% is more than enough given the distinct options.
    if math.isclose(theoretical_coefficient, answer_coefficient, rel_tol=0.01):
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        # Check for common mistakes that correspond to other options.
        
        # Mistake 1: Calculating the position of the first minimum (θ₁).
        # Coefficient would be z₁ / (2π) ≈ 0.610 (Option B)
        first_min_coeff = z1 / (2 * math.pi)
        
        # Mistake 2: Calculating the angular diameter of the central maximum (2 * θ₁).
        # Coefficient would be 2 * z₁ / (2π) = z₁ / π ≈ 1.220 (Option D)
        central_max_diam_coeff = z1 / math.pi

        return (f"Incorrect. The provided answer corresponds to a coefficient of {answer_coefficient}. "
                f"The correct theoretical coefficient for the angular distance between the first two minima is (z₂ - z₁) / (2π), "
                f"which calculates to approximately {theoretical_coefficient:.5f}. "
                f"This value matches option C (0.506). The provided answer 'C' is therefore correct, but if the reasoning led to a different option, that reasoning was flawed.\n"
                f"Common mistakes include:\n"
                f"- Calculating the position of the first minimum: coefficient ≈ {first_min_coeff:.5f} (matches option B).\n"
                f"- Calculating the angular diameter of the central maximum: coefficient ≈ {central_max_diam_coeff:.5f} (matches option D).")

# Run the check and print the result.
result = check_diffraction_answer()
print(result)