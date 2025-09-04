import math

def check_decay_probability():
    """
    Checks the correctness of the answer to the radioactive decay problem.

    The problem states:
    - Decay probability in 100 minutes is 32%.
    - Atom has survived for 50 minutes.
    - Question: What is the probability it will decay in the next 150 minutes?

    The solution relies on the memoryless property of exponential decay.
    """

    # --- Step 1: Define problem parameters ---
    p_decay_100_min = 0.32
    t1 = 100  # minutes
    t2 = 150  # minutes for the new probability calculation

    # The provided answer to check
    llm_answer_choice = 'B'
    options = {'A': 0.48, 'B': 0.44, 'C': 0.52, 'D': 0.40}

    # --- Step 2: Calculate the survival probability for the initial period ---
    # P(decay) = 1 - P(survival)
    # So, P(survival) = 1 - P(decay)
    p_survival_100_min = 1 - p_decay_100_min
    # This gives us S(100) = 0.68, where S(t) = exp(-lambda * t)

    # --- Step 3: Calculate the required probability using the memoryless property ---
    # The 50 minutes of survival are irrelevant. We need to find the probability
    # of decay in a 150-minute interval.
    # P(decay in 150 min) = 1 - S(150)
    # We know S(t) = exp(-lambda * t), so S(150) = exp(-lambda * 150)
    # We can relate S(150) to S(100):
    # S(150) = exp(-lambda * 100 * 1.5) = (exp(-lambda * 100))^1.5
    # S(150) = S(100) ^ (150/100)
    
    p_survival_150_min = p_survival_100_min ** (t2 / t1)
    
    # The probability of decay in 150 minutes is 1 minus the survival probability
    calculated_p_decay_150_min = 1 - p_survival_150_min

    # --- Step 4: Compare the calculated result with the given options ---
    # Find which option is numerically closest to our calculated result.
    closest_option = None
    min_difference = float('inf')

    for option_key, option_value in options.items():
        difference = abs(calculated_p_decay_150_min - option_value)
        if difference < min_difference:
            min_difference = difference
            closest_option = option_key
            
    # --- Step 5: Validate the LLM's answer ---
    # Check if the LLM's chosen option matches the closest option we found.
    if llm_answer_choice == closest_option:
        return "Correct"
    else:
        return (f"Incorrect. The reasoning provided in the final answer is correct, but the final choice is not.\n"
                f"Calculation Steps:\n"
                f"1. Probability of survival for 100 min = 1 - 0.32 = 0.68.\n"
                f"2. Due to the memoryless property, we calculate the probability of decay in a 150-minute interval.\n"
                f"3. Probability of survival for 150 min = (0.68)^(150/100) = 0.68^1.5 ≈ {p_survival_150_min:.4f}.\n"
                f"4. Probability of decay in 150 min = 1 - {p_survival_150_min:.4f} ≈ {calculated_p_decay_150_min:.4f} (or {calculated_p_decay_150_min*100:.2f}%).\n"
                f"5. The closest option to {calculated_p_decay_150_min*100:.2f}% is {closest_option} ({options[closest_option]*100}%), not {llm_answer_choice} ({options[llm_answer_choice]*100}%).")

# Run the check
result = check_decay_probability()
print(result)