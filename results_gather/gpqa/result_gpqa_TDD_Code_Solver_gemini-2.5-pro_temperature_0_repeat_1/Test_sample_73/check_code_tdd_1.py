import math

def check_answer():
    """
    This function verifies the solution to the radioactive decay problem.
    
    The problem states:
    1. Decay probability is 32% in 100 minutes.
    2. 50 minutes have passed without decay.
    3. What is the probability of decay in the next 150 minutes?

    The provided answer is A) 44%.
    """

    # --- Step 1: Define problem parameters and the given answer ---
    
    # Given information
    t1 = 100  # minutes
    p_decay_t1 = 0.32
    
    # Time interval for the question
    t2 = 150  # minutes
    
    # The provided answer is 'A', which corresponds to 44%
    llm_answer_choice = 'A'
    options = {'A': 0.44, 'B': 0.48, 'C': 0.52, 'D': 0.40}
    
    # --- Step 2: Calculate the decay constant (lambda) ---
    
    # The probability of an atom *not* decaying (surviving) is P_survival = exp(-lambda * t)
    # The probability of it decaying is P_decay = 1 - P_survival = 1 - exp(-lambda * t)
    # From the given info: 0.32 = 1 - exp(-lambda * 100)
    # Rearranging to solve for lambda:
    # exp(-lambda * 100) = 1 - 0.32 = 0.68
    # -lambda * 100 = ln(0.68)
    # lambda = -ln(0.68) / 100
    
    try:
        decay_constant_lambda = -math.log(1 - p_decay_t1) / t1
    except (ValueError, ZeroDivisionError) as e:
        return f"Error during calculation of the decay constant: {e}"

    # --- Step 3: Calculate the probability for the new time interval ---
    
    # Radioactive decay is a memoryless process. The 50 minutes that have passed
    # do not affect the future probability. The clock "resets".
    # We calculate the probability of decay over the next 150 minutes.
    # P_decay_t2 = 1 - exp(-lambda * t2)
    
    p_decay_t2 = 1 - math.exp(-decay_constant_lambda * t2)

    # --- Step 4: Verify the correctness of the LLM's answer ---
    
    # Find which option is numerically closest to our calculated probability.
    calculated_prob = p_decay_t2
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_prob))

    # Check if the LLM's chosen option is the one closest to the correct answer.
    if closest_option == llm_answer_choice:
        # The calculated probability (approx 43.93%) rounds to 44%, which is option A.
        # The LLM correctly chose A.
        return "Correct"
    else:
        # This block would execute if the LLM's answer was incorrect.
        reason = (
            f"Incorrect. The provided answer is {llm_answer_choice} ({options[llm_answer_choice]*100}%).\n"
            f"The calculation shows the probability of decay in 150 minutes is {calculated_prob:.4f}, "
            f"which is {calculated_prob*100:.2f}%.\n"
            f"This value is closest to option {closest_option} ({options[closest_option]*100}%)."
        )
        return reason

# Run the check
result = check_answer()
print(result)