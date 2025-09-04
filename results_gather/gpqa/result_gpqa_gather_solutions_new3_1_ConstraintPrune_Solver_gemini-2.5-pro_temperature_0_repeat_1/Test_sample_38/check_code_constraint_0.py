import math
from scipy import integrate

def check_answer():
    """
    Checks the correctness of the final answer by verifying the normalization condition.
    """
    # The options provided in the question
    options = {
        'A': 0.6,
        'B': 1.1,
        'C': 0.85,
        'D': 0.35
    }

    # The final answer provided by the LLM
    llm_answer_choice = 'C'
    
    if llm_answer_choice not in options:
        return f"Incorrect. The answer choice '{llm_answer_choice}' is not a valid option."

    # Get the numerical value of 'a' from the chosen answer
    a_val = options[llm_answer_choice]

    # Define the probability density function |ψ(x)|²
    # This is the function to be integrated.
    def probability_density(x, a):
        return (a**2) / (1 + x) + 0.25

    # Numerically calculate the integral from x=1 to x=3
    # The result should be 1 if the wave function is normalized.
    # integrate.quad returns a tuple (result, estimated_error)
    total_probability, _ = integrate.quad(probability_density, 1, 3, args=(a_val,))

    # The exact analytical solution for 'a' is sqrt(0.5 / ln(2))
    correct_a = math.sqrt(0.5 / math.log(2))

    # Check if the calculated total probability is close to 1.
    # We use a tolerance because the option a=0.85 is a rounded value.
    # The calculated probability for a=0.85 is ~1.0008.
    if math.isclose(total_probability, 1.0, rel_tol=1e-2):
        return "Correct"
    else:
        # Provide a detailed reason for the failure
        reason = (
            f"Incorrect. The normalization condition is not satisfied for the given answer.\n"
            f"The chosen answer is {llm_answer_choice}, which corresponds to a = {a_val}.\n"
            f"The normalization integral ∫[from 1 to 3] |ψ(x)|² dx must equal 1.\n"
            f"For a = {a_val}, the integral evaluates to approximately {total_probability:.4f}, which is not 1.\n"
            f"The correct value for 'a' is sqrt(0.5 / ln(2)) ≈ {correct_a:.4f}, making 0.85 the closest choice, but the final answer's reasoning should lead to this choice."
        )
        # Since the question is multiple choice and 0.85 is the closest option,
        # the choice 'C' is technically correct. The code checks the physical principle.
        # Let's re-evaluate the logic. The LLM correctly identified the closest option.
        # The check should confirm this.
        
        # Let's check which option gives the result closest to 1.
        results = {}
        for option, val in options.items():
            prob, _ = integrate.quad(probability_density, 1, 3, args=(val,))
            results[option] = prob
        
        closest_option = min(results, key=lambda k: abs(results[k] - 1))

        if llm_answer_choice == closest_option:
             return "Correct"
        else:
             return (f"Incorrect. The chosen option '{llm_answer_choice}' (a={a_val}) results in a total probability of {results[llm_answer_choice]:.4f}. "
                     f"The option '{closest_option}' (a={options[closest_option]}) is a better fit, with a total probability of {results[closest_option]:.4f}, which is closer to 1.")


# Run the check
print(check_answer())