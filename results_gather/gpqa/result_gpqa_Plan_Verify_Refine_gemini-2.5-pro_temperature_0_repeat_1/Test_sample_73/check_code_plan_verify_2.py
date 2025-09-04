import math

def check_radioactive_decay_answer():
    """
    This function verifies the answer to the radioactive decay problem.
    It calculates the correct probability and checks if the provided answer 'A' (44%)
    is the closest option among the choices.
    """
    # --- Problem Parameters ---
    # Given: 32% decay probability in 100 minutes.
    p_decay_100 = 0.32
    t1 = 100.0  # minutes

    # Question: What is the probability of decay in the next 150 minutes?
    # The 50 minutes that have passed are irrelevant due to the memoryless property.
    t2 = 150.0  # minutes

    # --- Multiple Choice Options ---
    options = {'A': 0.44, 'B': 0.48, 'C': 0.52, 'D': 0.40}
    chosen_answer_label = 'A'
    
    # --- Calculation ---
    # 1. Probability of survival for 100 minutes.
    p_survival_100 = 1 - p_decay_100

    # 2. The survival probability function is P_survival(t) = exp(-lambda * t).
    # This implies P_survival(t_b) = P_survival(t_a)^(t_b / t_a).
    # We can calculate the survival probability for 150 minutes without finding lambda explicitly.
    p_survival_150 = p_survival_100 ** (t2 / t1)

    # 3. The probability of decay in 150 minutes is 1 minus the survival probability.
    calculated_p_decay_150 = 1 - p_survival_150

    # --- Verification ---
    # Find which of the given options is closest to our calculated result.
    
    # Calculate the absolute difference between the calculated value and each option.
    differences = {label: abs(calculated_p_decay_150 - value) for label, value in options.items()}

    # Find the option label with the smallest difference.
    closest_option_label = min(differences, key=differences.get)

    # Check if the chosen answer ('A') is indeed the closest one.
    if closest_option_label == chosen_answer_label:
        return "Correct"
    else:
        # This part would execute if the answer was wrong.
        reason = (
            f"The provided answer '{chosen_answer_label}' ({options[chosen_answer_label]*100}%) is incorrect.\n"
            f"The calculated probability of decay in 150 minutes is {calculated_p_decay_150:.4f}, which is {calculated_p_decay_150*100:.2f}%.\n"
            f"The closest option to this value is '{closest_option_label}' ({options[closest_option_label]*100}%)."
        )
        return reason

# Run the check
result = check_radioactive_decay_answer()
print(result)