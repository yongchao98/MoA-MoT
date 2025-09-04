import math

def check_answer():
    """
    Checks the correctness of the answer to the radioactive decay problem.
    """
    # --- Problem Constraints & Data ---
    # Decay probability in 100 minutes is 32%
    p_decay_100 = 0.32
    t1 = 100  # minutes
    
    # We need to find the decay probability in the next 150 minutes.
    # Due to the memoryless property of exponential decay, the 50 minutes that have passed are irrelevant.
    t2 = 150  # minutes
    
    # The multiple-choice options provided in the question
    options = {'A': 0.40, 'B': 0.48, 'C': 0.44, 'D': 0.52}
    
    # The final answer given by the LLM to be checked
    llm_answer_choice = 'C'

    # --- Step 1: Calculate the survival probability for the known interval (100 min) ---
    # S(t) = 1 - P(t)
    s_100 = 1 - p_decay_100
    
    # --- Step 2: Calculate the survival probability for the target interval (150 min) ---
    # The survival function is S(t) = exp(-lambda*t).
    # Therefore, S(t2) = exp(-lambda*t2) = exp(-lambda*t1 * (t2/t1)) = (exp(-lambda*t1))^(t2/t1)
    # S(150) = S(100) ^ (150/100)
    s_150 = s_100 ** (t2 / t1)
    
    # --- Step 3: Calculate the decay probability for the target interval (150 min) ---
    p_decay_150 = 1 - s_150
    
    # --- Step 4: Find the closest option to the calculated probability ---
    closest_option = None
    min_difference = float('inf')
    
    for option, value in options.items():
        difference = abs(p_decay_150 - value)
        if difference < min_difference:
            min_difference = difference
            closest_option = option
            
    # --- Step 5: Verify the LLM's answer ---
    # Check if the reasoning and calculation lead to the provided answer.
    # The calculated value is ~43.93%, which is closest to 44% (Option C).
    if closest_option == llm_answer_choice:
        return "Correct"
    else:
        return (f"Incorrect. The calculated decay probability is {p_decay_150:.4f} ({p_decay_150:.2%}).\n"
                f"This value is closest to option {closest_option} ({options[closest_option]:.0%}).\n"
                f"The provided answer was {llm_answer_choice} ({options[llm_answer_choice]:.0%}), which is not the closest option.")

# Execute the check and print the result
result = check_answer()
print(result)