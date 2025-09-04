import math
import scipy.integrate

def check_correctness():
    """
    Checks the correctness of the final answer by:
    1. Analytically solving for 'a' and comparing it to the options.
    2. Numerically integrating the probability density with the given 'a' to see if it equals 1.
    """
    
    # The final answer from the LLM is D, which corresponds to a = 0.85
    final_answer_letter = 'D'
    options = {'A': 0.35, 'B': 1.1, 'C': 0.6, 'D': 0.85}
    
    if final_answer_letter not in options:
        return f"Incorrect: The final answer letter '{final_answer_letter}' is not a valid option."
        
    a_candidate = options[final_answer_letter]

    # --- Check 1: Analytical Solution ---
    # The normalization condition leads to the equation: a² * ln(2) + 0.5 = 1
    # Solving for 'a', we get: a = sqrt(0.5 / ln(2))
    try:
        a_exact = math.sqrt(0.5 / math.log(2))
    except Exception as e:
        return f"Incorrect: An error occurred during the analytical calculation: {e}"

    # Find which option is mathematically closest to the exact value
    closest_option_letter = min(options, key=lambda k: abs(options[k] - a_exact))

    if final_answer_letter != closest_option_letter:
        return (f"Incorrect: The analytical calculation gives a ≈ {a_exact:.4f}. "
                f"The closest option is {closest_option_letter} ({options[closest_option_letter]}), "
                f"but the provided answer was {final_answer_letter} ({a_candidate}).")

    # --- Check 2: Numerical Verification ---
    # Define the probability density function |ψ(x)|² = a²/(1+x) + 0.25
    def probability_density(x, a):
        return (a**2) / (1 + x) + 0.25

    # Integrate the function from x=1 to x=3 using the candidate 'a'
    try:
        # The args parameter passes the value of 'a' to the function
        integral_result, _ = scipy.integrate.quad(probability_density, 1, 3, args=(a_candidate,))
    except Exception as e:
        return f"Incorrect: An error occurred during the numerical integration check: {e}"

    # The result of the integral must be 1 for the normalization condition to be met.
    # We use a small tolerance because the option a=0.85 is a rounded value.
    if not math.isclose(integral_result, 1.0, rel_tol=1e-3):
        return (f"Incorrect: The normalization condition is not satisfied. "
                f"For the given answer a = {a_candidate}, the integral of the probability density "
                f"from 1 to 3 is {integral_result:.6f}, which is not equal to 1.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)