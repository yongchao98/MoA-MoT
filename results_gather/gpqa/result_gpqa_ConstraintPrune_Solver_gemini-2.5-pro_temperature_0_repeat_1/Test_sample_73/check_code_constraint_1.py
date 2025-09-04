import math

def check_answer():
    """
    This function checks the correctness of the LLM's answer to the radioactive decay problem.

    The problem involves two main steps:
    1. Calculate the decay constant (λ) from the initial condition (32% decay in 100 minutes).
    2. Use the memoryless property of exponential decay to calculate the probability of decay in a new interval (150 minutes).
    """

    # --- Define problem parameters ---
    # Initial condition: 32% probability of decay in 100 minutes.
    p_decay_initial = 0.32
    time_initial = 100  # minutes

    # Question: Probability of decay in the next 150 minutes.
    # Due to the memoryless property, the 50 minutes that have passed are irrelevant.
    time_future = 150  # minutes

    # --- LLM's Answer ---
    # The LLM's final answer is D, which corresponds to 44%.
    llm_answer_key = 'D'
    options = {'A': 0.52, 'B': 0.48, 'C': 0.40, 'D': 0.44}
    
    # --- Step 1: Calculate the decay constant (λ) ---
    # The probability of decay P(t) is given by P(t) = 1 - e^(-λt).
    # Rearranging to solve for λ:
    # e^(-λt) = 1 - P(t)
    # -λt = ln(1 - P(t))
    # λ = -ln(1 - P(t)) / t
    try:
        lambda_constant = -math.log(1 - p_decay_initial) / time_initial
    except (ValueError, TypeError) as e:
        return f"Error in calculation of the decay constant λ: {e}"

    # --- Step 2: Calculate the probability for the new time interval ---
    # P(decay in 150 min) = 1 - e^(-λ * 150)
    prob_decay_future = 1 - math.exp(-lambda_constant * time_future)

    # --- Step 3: Verify the LLM's answer ---
    # Find which option is numerically closest to our calculated probability.
    closest_option_key = min(options, key=lambda k: abs(options[k] - prob_decay_future))

    # Check if the LLM's chosen option is the correct one.
    if closest_option_key == llm_answer_key:
        # The logic and final answer choice are correct.
        # We can also check the intermediate calculation from the LLM's explanation (P ≈ 0.4393).
        llm_intermediate_value = 0.4393
        if abs(prob_decay_future - llm_intermediate_value) < 0.0001:
            return "Correct"
        else:
            # The intermediate value has a minor rounding difference, but the final conclusion is correct.
            # This is still considered correct.
            return (f"Correct. The final answer choice is correct. There is a minor rounding difference in the "
                    f"intermediate steps shown in the explanation. The code calculates the probability as "
                    f"{prob_decay_future:.6f}, while the explanation gives {llm_intermediate_value}. "
                    f"Both correctly round to 44%.")
    else:
        # The LLM's choice was incorrect.
        return (f"Incorrect. The calculated probability of decay in 150 minutes is {prob_decay_future:.4f} "
                f"({prob_decay_future*100:.2f}%). The closest option is {closest_option_key} "
                f"({options[closest_option_key]*100}%), but the provided answer was {llm_answer_key} "
                f"({options[llm_answer_key]*100}%).")

# Execute the check and print the result
result = check_answer()
# The provided answer is correct, but the code will return the more detailed "Correct" message
# because of minor rounding differences. For the final output, we simplify it to just "Correct".
if "Correct" in result:
    print("Correct")
else:
    print(result)