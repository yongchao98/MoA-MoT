import numpy as np
try:
    from scipy.special import jn_zeros
except ImportError:
    # Define a dummy function if scipy is not installed
    jn_zeros = None

def check_correctness():
    """
    Checks the correctness of the provided answer for the diffraction problem.
    It calculates the theoretical coefficient for the angular distance between
    the first two minima of an Airy pattern and compares it to the given answer.
    """
    # The provided answer is B, which corresponds to a coefficient of 0.506.
    answer_coefficient = 0.506

    # The first two non-zero roots of the Bessel function J1(x) are needed.
    # We can get them from scipy or use their well-known values.
    if jn_zeros:
        try:
            zeros = jn_zeros(1, 2)  # Get the first 2 positive roots of J1(x)
            z1, z2 = zeros[0], zeros[1]
        except Exception:
            # Fallback to hardcoded values if scipy fails for any reason
            z1 = 3.83170597
            z2 = 7.01558667
    else:
        # Hardcoded high-precision values if scipy is not available.
        z1 = 3.83170597
        z2 = 7.01558667

    # Calculate the theoretical coefficient for the angular distance: C = (z2 - z1) / (2 * pi)
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # Compare the calculated coefficient with the answer's coefficient.
    # A tolerance is used to account for rounding in the provided option.
    tolerance = 0.001
    if abs(calculated_coefficient - answer_coefficient) < tolerance:
        return "Correct"
    else:
        return (
            f"The answer is incorrect.\n"
            f"The problem describes diffraction from a circular aperture of radius 'a'.\n"
            f"The angular distance between the first two minima is given by Delta_theta = ((z2 - z1) / (2*pi)) * (lambda/a).\n"
            f"Using the first two zeros of the J1 Bessel function, z1 ≈ {z1:.4f} and z2 ≈ {z2:.4f}.\n"
            f"The calculated coefficient is ({z2:.4f} - {z1:.4f}) / (2 * pi) ≈ {calculated_coefficient:.4f}.\n"
            f"The coefficient from the selected answer (B) is {answer_coefficient}, but the calculated value is {calculated_coefficient:.4f}. "
            f"These values do not match."
        )

# Run the check
result = check_correctness()
print(result)