import math

def check_answer():
    """
    This function checks the correctness of the provided LLM answer for the radioactive decay problem.
    """
    # --- Problem Setup ---
    # Given: Probability of decay in 100 minutes is 32%
    p_decay_100 = 0.32
    
    # Question: Probability of decay in the next 150 minutes, given it has survived 50 minutes.
    # Options: A) 52%, B) 40%, C) 48%, D) 44%
    options = {'A': 52, 'B': 40, 'C': 48, 'D': 44}
    
    # LLM's final answer
    llm_answer_choice = 'D'

    # --- Step 1: Verify the core principle (Memoryless Property) ---
    # The LLM's reasoning correctly identifies that radioactive decay is a memoryless process.
    # This means the 50 minutes of survival are irrelevant.
    # The problem simplifies to finding the probability of decay in a 150-minute interval.
    # This reasoning is correct.

    # --- Step 2: Replicate the calculation ---
    # The probability of survival for 100 minutes is S(100) = 1 - P(100)
    s_100 = 1 - p_decay_100
    
    # The survival probability is S(t) = e^(-λt).
    # So, S(100) = e^(-100λ) = 0.68.
    
    # We need to find the decay probability for 150 minutes, P(150) = 1 - S(150).
    # S(150) = e^(-150λ) = (e^(-100λ))^(150/100) = (S(100))^1.5
    s_150 = s_100 ** 1.5
    
    # The decay probability for 150 minutes is:
    p_decay_150 = 1 - s_150
    
    # Convert to percentage
    calculated_percentage = p_decay_150 * 100

    # --- Step 3: Compare the calculation with the LLM's answer ---
    # The LLM's chosen option is D, which corresponds to 44%.
    llm_answer_value = options.get(llm_answer_choice)
    if llm_answer_value is None:
        return f"Incorrect. The LLM's answer choice '{llm_answer_choice}' is not a valid option."

    # Check if the LLM's chosen option is the closest to the calculated result.
    # This accounts for potential rounding in the problem statement or options.
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_percentage))

    if llm_answer_choice != closest_option:
        return (f"Incorrect. The calculated probability is approximately {calculated_percentage:.2f}%. "
                f"The closest option is {closest_option} ({options[closest_option]}%), but the LLM chose {llm_answer_choice} ({llm_answer_value}%).")

    # The LLM chose the correct option. Let's also check if its own calculation was accurate.
    # The LLM calculated ~43.93%, which rounds to 44%. Our calculation is ~43.93%. This is consistent.
    
    # The reasoning is sound, the calculation is correct, and the chosen option matches the calculation.
    return "Correct"

# Run the check
result = check_answer()
print(result)