import math

def check_answer():
    """
    This function verifies the solution to the atomic decay problem.
    
    The problem states:
    - Decay probability is 32% within 100 minutes.
    - 50 minutes have passed without decay.
    - Question: What is the probability of decay in the next 150 minutes?

    The core principle is that radioactive decay is a memoryless Poisson process.
    The probability of decay in a future interval is independent of how long the atom has already survived.
    
    The probability of decay P(t) in time t is given by P(t) = 1 - exp(-λt),
    where λ is the decay constant.
    """
    
    # --- Step 1: Define problem parameters ---
    p_initial = 0.32  # 32% decay probability
    t_initial = 100   # in 100 minutes
    t_future = 150    # we want the probability for the next 150 minutes
    
    # The 50 minutes that have passed are irrelevant due to the memoryless property of decay.

    # --- Step 2: Calculate the decay constant (λ) from the given data ---
    # p_initial = 1 - exp(-λ * t_initial)
    # exp(-λ * t_initial) = 1 - p_initial
    # -λ * t_initial = ln(1 - p_initial)
    # λ = -ln(1 - p_initial) / t_initial
    try:
        lambda_const = -math.log(1 - p_initial) / t_initial
    except ValueError:
        return "Error: Cannot calculate log of a non-positive number. Initial probability must be less than 1."

    # --- Step 3: Calculate the probability of decay in the future time interval ---
    # P(decay in t_future) = 1 - exp(-λ * t_future)
    prob_decay_future = 1 - math.exp(-lambda_const * t_future)
    
    # --- Step 4: Compare the calculated result with the provided answer ---
    # The provided answer is A) 44%
    llm_answer_value = 0.44
    
    # The options are A) 44%, B) 48%, C) 52%, D) 40%.
    # We check if our calculated result is closest to the chosen option A.
    # The difference between options is 4%, so a tolerance of 1% (0.01) is reasonable.
    if math.isclose(prob_decay_future, llm_answer_value, abs_tol=0.01):
        return "Correct"
    else:
        # If the answer is not correct, explain why.
        reason = (f"The calculated probability is {prob_decay_future:.4f}, which is {prob_decay_future*100:.2f}%. "
                  f"The provided answer is {llm_answer_value*100}%. These values are not close enough.\n"
                  f"Calculation steps:\n"
                  f"1. The decay constant λ is calculated from P(100 min) = 0.32. "
                  f"λ = -ln(1 - 0.32) / 100 ≈ {lambda_const:.6f}.\n"
                  f"2. The probability of decay in the next 150 minutes is P(150 min) = 1 - exp(-λ * 150).\n"
                  f"3. P(150 min) = 1 - exp(-{lambda_const:.6f} * 150) ≈ {prob_decay_future:.4f} or {prob_decay_future*100:.2f}%.")
        return reason

# Execute the check and print the result
result = check_answer()
print(result)