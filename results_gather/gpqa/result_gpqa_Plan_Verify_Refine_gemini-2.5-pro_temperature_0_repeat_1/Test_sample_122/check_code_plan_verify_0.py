import math

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the special relativity problem.
    """
    # --- Problem Constraints and Given Values ---
    v = 60000  # km/s
    dt0 = 50   # s (proper time in the ejecta's reference frame)
    c = 300000 # km/s (approximated speed of light)

    # --- LLM's Answer ---
    llm_answer_choice = 'A'
    options = {
        'A': 3060000,
        'B': 2940000,
        'C': 3000000,
        'D': 2880000
    }
    llm_answer_value = options.get(llm_answer_choice)

    if llm_answer_value is None:
        return f"Invalid answer choice '{llm_answer_choice}' provided."

    # --- Calculation ---
    # 1. Check for non-relativistic case (classical physics)
    classical_distance = v * dt0
    if classical_distance == llm_answer_value:
        return f"Incorrect. The answer corresponds to a non-relativistic calculation ({classical_distance} km), which is option C. The problem involves high speeds (20% the speed of light), so relativistic effects (time dilation) must be considered."

    # 2. Relativistic Calculation
    try:
        # Calculate beta (v/c)
        beta = v / c
        # Calculate Lorentz factor (gamma)
        gamma = 1 / math.sqrt(1 - beta**2)
        # Calculate time in the Galaxy's reference frame (dt)
        dt = gamma * dt0
        # Calculate distance in the Galaxy's reference frame
        calculated_distance = v * dt
    except (ValueError, ZeroDivisionError) as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # Find the option that is numerically closest to the calculated distance
    closest_option = min(options.keys(), key=lambda option: abs(options[option] - calculated_distance))

    # Check if the LLM's answer is the closest one
    if llm_answer_choice == closest_option:
        # Additionally, check if the relative error is small, confirming the option is a good approximation
        relative_error = abs(calculated_distance - llm_answer_value) / llm_answer_value
        if relative_error < 0.01: # 1% tolerance for rounding in the options
            return "Correct"
        else:
            # This case is unlikely for well-posed problems but is a good sanity check
            return f"Correct. The LLM chose the closest option ({llm_answer_choice}). However, the relative error between the calculated value ({calculated_distance:.2f} km) and the option value ({llm_answer_value} km) is {relative_error:.2%}, which is larger than expected."
    else:
        return f"Incorrect. The calculated relativistic distance is approximately {calculated_distance:.2f} km. The closest option is {closest_option} ({options[closest_option]} km), but the provided answer was {llm_answer_choice} ({llm_answer_value} km)."

# Run the check
result = check_correctness()
print(result)