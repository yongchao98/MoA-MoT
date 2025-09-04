import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the radioactive decay problem.
    """
    # Given information from the question
    # Probability of decay within 100 minutes is 32%
    p_decay_100 = 0.32
    time_1 = 100  # minutes
    
    # The question asks for the probability of decay in the next 150 minutes.
    # Due to the memoryless property of exponential decay, the 50 minutes that have passed are irrelevant.
    # We need to calculate the probability of decay in a 150-minute interval.
    time_2 = 150  # minutes

    # Step 1: Calculate the survival probability for 100 minutes.
    # S(t) = 1 - P(t)
    p_survival_100 = 1 - p_decay_100

    # Step 2: Calculate the survival probability for 150 minutes.
    # The relationship is S(t2) = S(t1)^(t2/t1).
    p_survival_150 = p_survival_100 ** (time_2 / time_1)

    # Step 3: Calculate the decay probability for 150 minutes.
    p_decay_150 = 1 - p_survival_150
    
    # Convert the result to a percentage.
    calculated_percentage = p_decay_150 * 100

    # The options provided in the question.
    options = {'A': 40, 'B': 52, 'C': 48, 'D': 44}
    
    # The answer provided by the LLM.
    llm_answer_letter = 'D'

    # Find which option is numerically closest to the calculated result.
    closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_percentage))

    # Check if the LLM's answer matches the closest calculated option.
    if llm_answer_letter == closest_option_letter:
        return "Correct"
    else:
        llm_answer_value = options.get(llm_answer_letter)
        closest_option_value = options[closest_option_letter]
        return (f"Incorrect. The calculated probability is approximately {calculated_percentage:.2f}%. "
                f"The closest option is {closest_option_value}% ({closest_option_letter}). "
                f"The provided answer chose {llm_answer_value}% ({llm_answer_letter}).")

# Run the check
result = check_correctness()
print(result)