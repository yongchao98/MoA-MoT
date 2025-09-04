import numpy as np
try:
    from scipy.special import jn_zeros
except ImportError:
    # Provide a dummy function if scipy is not available, so the code can still run and give a meaningful error.
    def jn_zeros(n, nt):
        raise ImportError("The 'scipy' library is required to run this check.")

def check_correctness():
    """
    This function checks the correctness of the provided answer to the physics problem.

    The problem asks for the angular distance between the first two diffraction minima
    for a circular aperture of radius 'a'.

    The physics behind the solution is as follows:
    1. An N-sided polygon with N->infinity and constant apothem 'a' becomes a circle of radius 'a'.
    2. The diffraction pattern is an Airy pattern.
    3. The angular position of the m-th minimum (for small angles) is given by:
       θ_m = (z_m * λ) / (2 * π * a)
       where z_m is the m-th non-trivial zero of the first-order Bessel function J₁(x).
    4. The angular distance between the first two minima is Δθ = θ₂ - θ₁.
       Δθ = [(z₂ - z₁) / (2 * π)] * (λ / a)
    5. The code calculates the coefficient C = (z₂ - z₁) / (2 * π) and compares it
       to the coefficient from the selected answer.
    """
    try:
        # Get the first two non-trivial zeros of the Bessel function J1(x).
        # jn_zeros(n, nt) computes the first nt zeros of J_n(x).
        # n=1 for J1, nt=2 for the first two zeros.
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]
        z2 = zeros[1]
    except ImportError:
        return "Could not verify the answer because the 'scipy' library is not installed. Please install it using 'pip install scipy'."
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Calculate the correct coefficient and the coefficients for common distractors ---

    # The correct coefficient for the angular distance between the first two minima
    correct_coefficient = (z2 - z1) / (2 * np.pi)

    # Coefficient for the position of the first minimum (a common distractor, option A)
    # θ₁ = [z₁ / (2π)] * (λ/a)
    distractor_coeff_1 = z1 / (2 * np.pi)

    # Coefficient for the angular diameter of the central bright disk (another common distractor, option B)
    # 2 * θ₁ = [2 * z₁ / (2π)] * (λ/a)
    distractor_coeff_2 = 2 * z1 / (2 * np.pi)

    # --- Check the provided answer ---

    # The question's options and the LLM's final answer
    # A) 0.610 λ / a
    # B) 1.220 λ /a
    # C) 0.500 λ / a
    # D) 0.506 λ / a
    # The final answer provided is <<<D>>>.
    llm_answer_coefficient = 0.506

    # Use a tolerance for floating-point comparisons. 1e-3 is reasonable for 3-decimal-place options.
    tolerance = 1e-3

    # 1. Check if the answer is correct
    if np.isclose(llm_answer_coefficient, correct_coefficient, atol=tolerance):
        return "Correct"

    # 2. If not correct, check if it matches a known distractor
    if np.isclose(llm_answer_coefficient, distractor_coeff_1, atol=tolerance):
        return (f"Incorrect. The selected answer's coefficient ({llm_answer_coefficient}) corresponds to the angular position of the *first* minimum (calculated as {distractor_coeff_1:.4f}), "
                f"not the angular distance between the first two minima.")

    if np.isclose(llm_answer_coefficient, distractor_coeff_2, atol=tolerance):
        return (f"Incorrect. The selected answer's coefficient ({llm_answer_coefficient}) corresponds to the angular *diameter* of the central bright disk (calculated as {distractor_coeff_2:.4f}), "
                f"not the angular distance between the first two minima.")

    # 3. If it's neither correct nor a known distractor, it's simply numerically wrong.
    return (f"Incorrect. The selected answer's coefficient ({llm_answer_coefficient}) does not match the calculated correct coefficient for the angular distance "
            f"between the first two minima, which is approximately {correct_coefficient:.4f}.")

# Execute the check and print the result
print(check_correctness())