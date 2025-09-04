import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer to the radioactive decay problem.
    """
    # --- Problem Definition ---
    # Given: Probability of decay in 100 minutes is 32%.
    p_decay_t1 = 0.32
    t1 = 100.0  # minutes

    # Question: What is the probability of decay in the next 150 minutes?
    # The 50 minutes that have passed are irrelevant due to the memoryless property of radioactive decay.
    t2 = 150.0  # minutes

    # --- Calculation ---
    # Step 1: Find the survival probability for the first time interval (t1).
    # The survival probability S(t) = 1 - P(t), where P(t) is the decay probability.
    # S(t) is also given by e^(-lambda*t).
    s_t1 = 1 - p_decay_t1  # This is e^(-lambda * 100)

    # Step 2: Calculate the survival probability for the second time interval (t2).
    # We need to find S(t2) = e^(-lambda * t2).
    # We can write S(t2) = e^(-lambda * t1 * (t2/t1)) = (e^(-lambda * t1))^(t2/t1).
    # So, S(150) = S(100)^(150/100).
    s_t2 = s_t1 ** (t2 / t1)

    # Step 3: Calculate the decay probability for the second time interval (t2).
    calculated_p_decay_t2 = 1 - s_t2

    # --- Verification ---
    # The options provided in the question.
    options = {
        'A': 0.40,
        'B': 0.52,
        'C': 0.44,
        'D': 0.48
    }

    # The final answer provided by the LLM.
    llm_answer_key = 'C'

    # Find which option is numerically closest to our calculated result.
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_p_decay_t2))

    # Check if the LLM's answer key matches the key of the closest option.
    if llm_answer_key == closest_option_key:
        # The answer is correct as it points to the right option.
        # We can also check if the calculation is reasonably close to the option's value.
        if math.isclose(calculated_p_decay_t2, options[llm_answer_key], rel_tol=0.01):
            return "Correct"
        else:
            # This case is for when the letter is correct but the value is far off, which is unlikely here.
            return (f"Incorrect. The calculated probability is {calculated_p_decay_t2:.4f} ({calculated_p_decay_t2*100:.2f}%). "
                    f"While this is closest to option {closest_option_key} ({options[closest_option_key]*100:.0f}%), "
                    f"the value is not close enough to the provided answer's value.")
    else:
        return (f"Incorrect. The calculated probability of decay in 150 minutes is {calculated_p_decay_t2:.4f} "
                f"({calculated_p_decay_t2*100:.2f}%). This value is closest to option {closest_option_key} "
                f"({options[closest_option_key]*100:.0f}%). The provided answer was {llm_answer_key} "
                f"({options[llm_answer_key]*100:.0f}%).")

# Run the check and print the result.
print(check_correctness())