import numpy as np
from scipy.integrate import quad

def check_answer():
    """
    Checks the correctness of the answer by solving for 'a' analytically and
    verifying it against the provided options and the normalization constraint.
    """
    # --- Step 1: Solve for 'a' analytically ---
    # The normalization condition is: ∫[1 to 3] (a^2 / (1 + x) + 0.25) dx = 1
    # Let's solve the integral:
    # ∫(a^2/(1+x) + 0.25) dx = a^2 * ln(1+x) + 0.25*x
    # Evaluating from 1 to 3:
    # [a^2*ln(1+3) + 0.25*3] - [a^2*ln(1+1) + 0.25*1]
    # = (a^2*ln(4) + 0.75) - (a^2*ln(2) + 0.25)
    # = a^2 * (ln(4) - ln(2)) + 0.5
    # = a^2 * ln(4/2) + 0.5
    # = a^2 * ln(2) + 0.5
    #
    # Now, set the result to 1:
    # a^2 * ln(2) + 0.5 = 1
    # a^2 * ln(2) = 0.5
    # a^2 = 0.5 / ln(2)
    # a = sqrt(0.5 / ln(2))

    try:
        a_analytical = np.sqrt(0.5 / np.log(2))
    except Exception as e:
        return f"An error occurred during analytical calculation: {e}"

    # --- Step 2: Check the provided answer <<<D>>> ---
    # The options are A) 1.1, B) 0.6, C) 0.35, D) 0.85
    # The provided answer is D, which corresponds to a = 0.85.
    answer_label = 'D'
    answer_a_value = 0.85

    # Check if the provided answer's 'a' value is the closest to the analytical solution.
    candidates = {
        "A": 1.1,
        "B": 0.6,
        "C": 0.35,
        "D": 0.85
    }
    
    # Find the candidate with the minimum absolute difference from the analytical value.
    best_label = min(candidates, key=lambda k: abs(candidates[k] - a_analytical))

    if best_label != answer_label:
        return (f"Incorrect. The analytical value for 'a' is approximately {a_analytical:.4f}. "
                f"The closest option is {best_label} (a={candidates[best_label]}), "
                f"but the provided answer was {answer_label} (a={answer_a_value}).")

    # --- Step 3: Verify the normalization constraint for the given answer ---
    # This confirms that the chosen option correctly normalizes the wave function.
    def probability_density(x, a):
        return (a**2 / (1 + x)) + 0.25

    # Perform numerical integration for the given answer's 'a' value.
    total_probability, _ = quad(probability_density, 1, 3, args=(answer_a_value,))

    # The total probability must be approximately 1. We use a tolerance
    # to account for the multiple-choice nature of the question.
    if not np.isclose(total_probability, 1.0, atol=1e-2):
        return (f"Incorrect. The normalization constraint is not satisfied for a = {answer_a_value}. "
                f"The integral of the probability density from x=1 to x=3 is {total_probability:.4f}, which is not close to 1.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)