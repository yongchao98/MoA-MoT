import math

def check_correctness():
    """
    Checks the correctness of the provided answer for the radioactive decay problem.

    The problem states:
    - Decay probability in 100 minutes is 32%.
    - Atom has survived for 50 minutes.
    - Question: What is the probability of decay in the next 150 minutes?

    The key principle is the "memoryless" property of exponential decay, which means
    the 50 minutes of survival are irrelevant. The problem simplifies to finding the
    decay probability over a 150-minute interval.
    """

    # Given parameters
    prob_decay_100_min = 0.32
    time1 = 100.0
    time2 = 150.0

    # The options provided in the question
    options = {'A': 0.40, 'B': 0.48, 'C': 0.44, 'D': 0.52}
    
    # The final answer provided by the LLM
    llm_answer_choice = 'C'

    # --- Step 1: Calculate the survival probability for 100 minutes ---
    # The probability of survival S(t) is 1 - P(decay, t)
    prob_survival_100_min = 1 - prob_decay_100_min
    # This should be 1 - 0.32 = 0.68

    # --- Step 2: Calculate the survival probability for 150 minutes ---
    # The survival probability is given by S(t) = exp(-lambda * t).
    # Therefore, S(100) = exp(-lambda * 100) = 0.68.
    # We need to find S(150) = exp(-lambda * 150).
    # We can write S(150) = exp(-lambda * 100 * 1.5) = (exp(-lambda * 100))^1.5
    # So, S(150) = S(100) ^ (150/100)
    prob_survival_150_min = prob_survival_100_min ** (time2 / time1)

    # --- Step 3: Calculate the decay probability for 150 minutes ---
    calculated_prob_decay_150_min = 1 - prob_survival_150_min

    # --- Step 4: Find the closest option to the calculated probability ---
    closest_option = None
    min_difference = float('inf')

    for option_key, option_value in options.items():
        difference = abs(calculated_prob_decay_150_min - option_value)
        if difference < min_difference:
            min_difference = difference
            closest_option = option_key

    # --- Step 5: Verify the LLM's answer ---
    # Check if the LLM's chosen option is indeed the closest one.
    if llm_answer_choice == closest_option:
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability of decay in 150 minutes is {calculated_prob_decay_150_min:.4f} "
                f"({calculated_prob_decay_150_min*100:.2f}%). "
                f"The closest option is {closest_option} ({options[closest_option]*100}%), "
                f"but the provided answer was {llm_answer_choice}.")

# Execute the check and print the result
result = check_correctness()
print(result)