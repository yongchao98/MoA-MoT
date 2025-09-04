import math

def check_decay_probability():
    """
    Checks the correctness of the given answer for the radioactive decay problem.
    """
    # --- Problem Parameters ---
    # Given: Probability of decay within 100 minutes is 32%
    p_decay_100 = 0.32
    t1 = 100.0  # minutes

    # Question: Given it has survived 50 minutes, what is the probability
    # it will decay in the next 150 minutes?
    t_future = 150.0

    # The answer to check from the LLM (Option A)
    llm_answer = 0.44  # 44%

    # --- Step 1: Calculate the decay constant (lambda) ---
    # The probability of survival up to time t is P_survival(t) = exp(-lambda * t).
    # The probability of decay within time t is P_decay(t) = 1 - P_survival(t) = 1 - exp(-lambda * t).
    # From the given info: 0.32 = 1 - exp(-lambda * 100)
    # This means: exp(-lambda * 100) = 1 - 0.32 = 0.68
    # Taking the natural log: -lambda * 100 = ln(0.68)
    # Solving for lambda:
    try:
        lambda_const = -math.log(1 - p_decay_100) / t1
    except (ValueError, TypeError):
        return "Error in initial data: Probability must be between 0 and 1."

    # --- Step 2: Calculate the required conditional probability ---
    # Radioactive decay is a memoryless process. The probability of decay in a future
    # interval is independent of how long the atom has already survived.
    # P(decay in next 150 min | survived 50 min) = P(decay in 150 min)
    # So, we calculate the probability of a fresh atom decaying in 150 minutes.
    
    calculated_prob = 1 - math.exp(-lambda_const * t_future)

    # --- Step 3: Compare the calculated result with the provided answer ---
    # Use a small tolerance for floating-point comparison.
    tolerance = 0.001  # 0.1% tolerance

    if abs(calculated_prob - llm_answer) < tolerance:
        return "Correct"
    else:
        reason = (
            f"The provided answer is {llm_answer*100:.0f}%, but the calculated probability is {calculated_prob*100:.2f}%.\n\n"
            f"Reasoning:\n"
            f"1. The decay constant λ is calculated from P(decay, 100 min) = 0.32. This gives λ = -ln(1 - 0.32) / 100 ≈ {lambda_const:.6f}.\n"
            f"2. Radioactive decay is a memoryless process. The fact that the atom has survived for 50 minutes is irrelevant to its future decay probability.\n"
            f"3. The question is equivalent to finding the probability that a new atom will decay within 150 minutes.\n"
            f"4. P(decay, 150 min) = 1 - exp(-λ * 150) = 1 - exp(-{lambda_const:.6f} * 150) ≈ {calculated_prob:.4f}.\n"
            f"5. This is a probability of {calculated_prob*100:.2f}%, which matches the 44% from option A, but the provided answer was likely a different option."
        )
        # This part of the code will only be reached if the LLM answer is NOT 0.44
        # Since the provided answer is A (44%), the code will return "Correct".
        # If the provided answer was, for example, B (48%), this block would execute.
        return f"Incorrect. The calculated probability is {calculated_prob*100:.2f}%, not {llm_answer*100:.0f}%."


# Execute the check
result = check_decay_probability()
print(result)