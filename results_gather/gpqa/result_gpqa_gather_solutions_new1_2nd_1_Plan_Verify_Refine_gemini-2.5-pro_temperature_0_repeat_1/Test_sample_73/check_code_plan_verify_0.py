import math

def check_decay_probability():
    """
    Checks the correctness of the final answer for the atomic decay problem.
    """
    # --- Problem Setup ---
    # Given: Decay probability in 100 minutes is 32%.
    p_decay_100 = 0.32
    
    # The atom has survived 50 minutes. We need the probability of decay in the next 150 minutes.
    # Due to the "memoryless" property of exponential decay, the 50 minutes are irrelevant.
    # The problem is to find the decay probability for a 150-minute interval.
    t_target = 150.0
    t_given = 100.0
    
    # The multiple-choice options provided in the question.
    options = {'A': 0.40, 'B': 0.44, 'C': 0.52, 'D': 0.48}
    
    # The final answer provided by the LLM.
    llm_final_answer = 'B'

    # --- Calculation ---
    # 1. Calculate the survival probability for 100 minutes.
    # P_survival(t) = 1 - P_decay(t)
    s_100 = 1 - p_decay_100
    
    # 2. Use the survival probability to find the probability for 150 minutes.
    # The survival function is S(t) = exp(-lambda * t).
    # Therefore, S(t_target) = exp(-lambda * t_target)
    #                        = exp(-lambda * t_given * (t_target / t_given))
    #                        = (exp(-lambda * t_given)) ^ (t_target / t_given)
    #                        = S(t_given) ^ (t_target / t_given)
    s_150 = s_100 ** (t_target / t_given)
    
    # 3. Calculate the decay probability for 150 minutes.
    p_decay_150 = 1 - s_150
    
    # --- Verification ---
    # Find which option is closest to our calculated result.
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - p_decay_150))
    
    # Check if the LLM's final answer matches the calculated closest option.
    if closest_option == llm_final_answer:
        # Further check if the calculation itself is sound.
        # The expected value is ~0.4393.
        if math.isclose(p_decay_150, 0.4393, abs_tol=0.0001):
            return "Correct"
        else:
            return f"The final answer {llm_final_answer} is correct by chance, but the underlying calculation is flawed. Calculated probability: {p_decay_150:.4f}"
    else:
        return (f"Incorrect. The calculated decay probability is {p_decay_150*100:.2f}%. "
                f"This is closest to option {closest_option} ({options[closest_option]*100}%), "
                f"but the provided answer was {llm_final_answer}.")

# Run the check and print the result.
result = check_decay_probability()
print(result)