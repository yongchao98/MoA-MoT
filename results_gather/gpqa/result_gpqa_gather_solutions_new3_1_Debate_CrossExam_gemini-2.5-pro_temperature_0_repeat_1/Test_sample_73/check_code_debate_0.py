import math

def check_decay_probability():
    """
    Checks the correctness of the answer to the radioactive decay problem.
    """
    # --- Problem Constraints & Given Data ---
    # Decay probability in 100 minutes is 32%
    p_decay_100 = 0.32
    t1 = 100  # minutes

    # We need to find the decay probability in the next 150 minutes.
    # Due to the memoryless property of radioactive decay, the 50 minutes
    # that have passed are irrelevant.
    t2 = 150  # minutes

    # The options provided in the question
    options = {'A': 0.52, 'B': 0.48, 'C': 0.44, 'D': 0.40}
    
    # The final answer from the LLM to be checked
    llm_answer_choice = 'C'

    # --- Step 1: Calculate the survival probability for 100 minutes ---
    # The survival probability S(t) is 1 - P(decay, t).
    # S(t) is also given by the formula e^(-lambda * t).
    s_100 = 1 - p_decay_100
    # s_100 = e^(-lambda * 100) = 0.68

    # --- Step 2: Calculate the survival probability for 150 minutes ---
    # We need to find S(150) = e^(-lambda * 150).
    # We can relate S(150) to S(100) without solving for lambda:
    # S(150) = e^(-lambda * 150) = e^(-lambda * 100 * 1.5) = (e^(-lambda * 100))^1.5
    s_150 = s_100 ** 1.5

    # --- Step 3: Calculate the decay probability for 150 minutes ---
    p_decay_150 = 1 - s_150

    # --- Step 4: Check the result against the provided answer ---
    # Find which option is numerically closest to our calculated probability.
    min_diff = float('inf')
    closest_option = None
    for option_key, option_value in options.items():
        diff = abs(p_decay_150 - option_value)
        if diff < min_diff:
            min_diff = diff
            closest_option = option_key

    # Verify if the LLM's chosen answer matches the closest option.
    if closest_option == llm_answer_choice:
        return "Correct"
    else:
        reason = (f"Incorrect. The calculated decay probability is {p_decay_150:.4f} ({p_decay_150*100:.2f}%).\n"
                  f"This value is closest to option {closest_option} ({options[closest_option]*100}%), "
                  f"but the provided answer was {llm_answer_choice} ({options[llm_answer_choice]*100}%).")
        return reason

# Run the check
result = check_decay_probability()
print(result)