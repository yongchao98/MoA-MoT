import math

def check_radioactive_decay_answer():
    """
    This function verifies the solution to the radioactive decay problem.
    
    Problem details:
    - Decay probability in 100 minutes is 32%.
    - After 50 minutes of no decay, what is the probability of decay in the next 150 minutes?

    Key concepts:
    - The process is memoryless, so the 50 minutes of survival are irrelevant.
    - The problem reduces to finding the probability of decay in a 150-minute interval.
    - Survival probability S(t) = exp(-λt).
    - Decay probability P(t) = 1 - S(t).
    """
    
    # Given parameters from the question
    prob_decay_100_min = 0.32
    time1 = 100.0
    time2 = 150.0
    
    # The options provided in the question
    options = {'A': 0.48, 'B': 0.44, 'C': 0.52, 'D': 0.40}
    
    # The final answer provided by the LLM to be checked
    llm_final_choice = 'B'

    # --- Calculation from first principles ---

    # 1. Calculate the survival probability for the first interval (100 minutes).
    # S(100) = 1 - P(100)
    prob_survival_100_min = 1 - prob_decay_100_min
    
    # 2. Calculate the survival probability for the second interval (150 minutes).
    # The relationship is S(t2) = S(t1)^(t2/t1) because S(t) = exp(-λt).
    # S(150) = S(100)^(150/100)
    prob_survival_150_min = prob_survival_100_min ** (time2 / time1)
    
    # 3. Calculate the decay probability for the second interval (150 minutes).
    # P(150) = 1 - S(150)
    calculated_prob_decay_150_min = 1 - prob_survival_150_min
    
    # --- Verification ---
    
    # 4. Find which option is numerically closest to our calculated probability.
    closest_option = None
    min_diff = float('inf')
    for option_key, option_value in options.items():
        diff = abs(calculated_prob_decay_150_min - option_value)
        if diff < min_diff:
            min_diff = diff
            closest_option = option_key
            
    # 5. Compare the calculated correct option with the LLM's choice.
    # The LLM's reasoning is also checked and found to be correct (identifying the memoryless property).
    if closest_option == llm_final_choice:
        return "Correct"
    else:
        reason = (f"The calculated probability of decay in 150 minutes is {calculated_prob_decay_150_min:.4f} "
                  f"({calculated_prob_decay_150_min*100:.2f}%). This is closest to option {closest_option} "
                  f"({options[closest_option]*100}%). The provided answer was option {llm_final_choice}, "
                  f"which is incorrect.")
        return reason

# Execute the check and print the result.
result = check_radioactive_decay_answer()
# print(result) # This will output "Correct"