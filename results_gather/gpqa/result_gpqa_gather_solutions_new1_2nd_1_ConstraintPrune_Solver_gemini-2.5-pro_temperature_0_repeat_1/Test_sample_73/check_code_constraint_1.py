import math

def check_decay_probability():
    """
    Checks the correctness of the provided answer for the atomic decay problem.

    The core logic relies on the "memoryless" property of exponential decay.
    The 50 minutes that have passed are irrelevant to the future probability.
    """
    # --- Problem Constraints and Data ---
    # Decay probability in 100 minutes is 32%
    p_decay_100 = 0.32
    t1 = 100  # minutes

    # We need to find the decay probability in the next 150 minutes
    t2 = 150  # minutes

    # --- Options as provided in the final answer block ---
    # A) 40%, B) 48%, C) 44%, D) 52%
    options = {'A': 40, 'B': 48, 'C': 44, 'D': 52}
    
    # The final answer provided to be checked
    provided_answer_letter = 'C'

    # --- Step-by-step Calculation ---

    # 1. The probability of decay P(t) is given by 1 - e^(-λt).
    #    The probability of survival S(t) is e^(-λt).
    #    Therefore, S(t) = 1 - P(t).
    s_100 = 1 - p_decay_100
    if not math.isclose(s_100, 0.68):
        return f"Calculation Error: Survival probability for 100 minutes should be 1 - 0.32 = 0.68, but got {s_100}."

    # 2. We need to find the decay probability for 150 minutes, P(150).
    #    P(150) = 1 - S(150)
    #    We can find S(150) from S(100) using the property of exponents:
    #    S(t2) = e^(-λ*t2) = e^(-λ*t1 * (t2/t1)) = (e^(-λ*t1))^(t2/t1) = S(t1)^(t2/t1)
    s_150 = s_100 ** (t2 / t1)

    # 3. Calculate the final decay probability for 150 minutes.
    p_decay_150 = 1 - s_150

    # Convert to a percentage for comparison with options.
    calculated_percentage = p_decay_150 * 100

    # --- Verification ---

    # Find the option closest to our calculated result.
    closest_option_letter = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_percentage - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_letter = letter

    # Check if the calculated closest option matches the provided answer.
    if closest_option_letter == provided_answer_letter:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{provided_answer_letter}' is incorrect.\n"
            f"The core principle is the 'memoryless' property of radioactive decay, so the initial 50 minutes are irrelevant.\n"
            f"1. The survival probability for 100 minutes is S(100) = 1 - 0.32 = 0.68.\n"
            f"2. The survival probability for 150 minutes is S(150) = S(100)^(150/100) = 0.68^1.5 ≈ {s_150:.4f}.\n"
            f"3. The decay probability for 150 minutes is P(150) = 1 - S(150) ≈ 1 - {s_150:.4f} = {p_decay_150:.4f}.\n"
            f"4. In percentage, the calculated probability is {calculated_percentage:.2f}%.\n"
            f"5. The closest option to {calculated_percentage:.2f}% is '{closest_option_letter}' ({options[closest_option_letter]}%).\n"
            f"6. The provided answer was '{provided_answer_letter}', but the calculated correct option is '{closest_option_letter}'."
        )
        return reason

# Execute the check and print the result
print(check_decay_probability())