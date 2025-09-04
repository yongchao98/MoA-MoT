import math

def check_answer():
    """
    Checks the correctness of the provided answer to the radioactive decay problem.
    """
    # --- Problem Parameters ---
    # Probability of decay within 100 minutes is 32%
    p_decay_100 = 0.32
    t1 = 100  # minutes

    # We need to find the probability of decay in the next 150 minutes.
    # Due to the memoryless property of radioactive decay, the 50 minutes that have passed are irrelevant.
    # The problem is equivalent to finding the probability of decay in a 150-minute interval.
    t2 = 150  # minutes

    # The provided answer is B, which corresponds to 44%
    options = {'A': 0.52, 'B': 0.44, 'C': 0.40, 'D': 0.48}
    llm_answer_char = 'B'
    
    # --- Calculation ---
    # The probability of survival for a time t is P_survival(t) = exp(-lambda * t)
    # The probability of decay is P_decay(t) = 1 - P_survival(t)

    # From the given information:
    # P_decay(100) = 1 - exp(-lambda * 100) = 0.32
    # Therefore, the survival probability for 100 minutes is:
    p_survival_100 = 1 - p_decay_100
    # p_survival_100 = exp(-lambda * 100) = 0.68

    # We need to find P_decay(150) = 1 - exp(-lambda * 150)
    # We can express exp(-lambda * 150) in terms of exp(-lambda * 100):
    # exp(-lambda * 150) = exp(-lambda * 100 * 1.5) = (exp(-lambda * 100))^1.5
    
    # So, the survival probability for 150 minutes is:
    p_survival_150 = p_survival_100 ** (t2 / t1)
    
    # And the decay probability for 150 minutes is:
    calculated_p_decay_150 = 1 - p_survival_150

    # --- Verification ---
    # Find which option is closest to the calculated probability
    closest_option = None
    min_diff = float('inf')

    for option_char, option_value in options.items():
        diff = abs(calculated_p_decay_150 - option_value)
        if diff < min_diff:
            min_diff = diff
            closest_option = option_char

    # Check if the LLM's answer matches the closest option
    if llm_answer_char == closest_option:
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability of decay in 150 minutes is {calculated_p_decay_150:.4f} (or {calculated_p_decay_150*100:.2f}%). "
                f"The closest option is {closest_option} ({options[closest_option]*100}%), but the provided answer was {llm_answer_char} ({options[llm_answer_char]*100}%).")

# Run the check
result = check_answer()
print(result)