import numpy as np

def check_decay_probability():
    """
    This function verifies the solution to the atomic decay problem.
    It calculates the exact probability based on the principles of exponential decay
    and compares it to the provided answer.
    """
    # --- 1. Define problem parameters from the question ---
    # Probability of decay within the first 100 minutes
    p_decay_initial = 0.32
    # Time interval for the initial probability
    t_initial = 100.0  # minutes
    # Time the atom has already survived
    t_survived = 50.0  # minutes
    # The future time interval in question
    t_future = 150.0  # minutes

    # The multiple-choice options provided in the question
    options = {'A': 0.44, 'B': 0.48, 'C': 0.52, 'D': 0.40}
    # The answer given by the LLM
    llm_answer_key = 'A'

    # --- 2. Calculate the decay constant (lambda) ---
    # The probability of decay P(t) in time t is given by P(t) = 1 - exp(-lambda * t).
    # We can rearrange this to find lambda:
    # P(t) = 1 - exp(-lambda * t)
    # 1 - P(t) = exp(-lambda * t)
    # ln(1 - P(t)) = -lambda * t
    # lambda = -ln(1 - P(t)) / t
    try:
        lambda_rate = -np.log(1 - p_decay_initial) / t_initial
    except ValueError:
        return "Error: Cannot calculate log of a non-positive number. The initial probability must be less than 1."

    # --- 3. Calculate the required probability using the memoryless property ---
    # Radioactive decay is a Poisson process, which is memoryless.
    # This means the probability of an event occurring in a future interval is independent
    # of how much time has already passed without the event.
    # Therefore, P(decay in next 150 min | survived 50 min) = P(decay in 150 min).
    # We use the same formula P(t) = 1 - exp(-lambda * t) with t = 150 minutes.
    
    calculated_prob = 1 - np.exp(-lambda_rate * t_future)

    # --- 4. Verify the LLM's answer ---
    # Check if the LLM's chosen option is the correct one.
    # First, find which option is numerically closest to our calculated result.
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_prob))

    # Second, check if the LLM's answer matches the closest option.
    if llm_answer_key != closest_option_key:
        return (f"Incorrect. The calculated probability is {calculated_prob:.4f} ({calculated_prob:.2%}). "
                f"This corresponds to option {closest_option_key} ({options[closest_option_key]:.2%}), "
                f"but the provided answer was {llm_answer_key}.")

    # Third, ensure the calculated value is reasonably close to the option's value.
    # A tolerance of 0.005 (0.5%) is reasonable for rounding.
    if not np.isclose(calculated_prob, options[llm_answer_key], atol=0.005):
        return (f"Incorrect. The calculated probability is {calculated_prob:.4f} ({calculated_prob:.2%}). "
                f"This is not sufficiently close to the value of the chosen answer {llm_answer_key} ({options[llm_answer_key]:.2%}).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_decay_probability()
print(result)