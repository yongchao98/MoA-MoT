import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer to the special relativity problem.
    """
    # --- Define problem parameters ---
    v = 60000  # Relative velocity in km/s
    dt0 = 50   # Proper time in seconds
    c = 300000 # Approximate speed of light in km/s

    # --- Define the options and the proposed answer ---
    options = {
        'A': 2940000,
        'B': 3000000,
        'C': 3060000,
        'D': 2880000
    }
    proposed_answer_letter = 'C'

    # --- Perform the physics calculation ---
    try:
        # Calculate the Lorentz factor (gamma)
        beta_squared = (v / c)**2
        gamma = 1 / math.sqrt(1 - beta_squared)

        # Calculate the time elapsed in the Galaxy's frame (dilated time)
        dt_galaxy = gamma * dt0

        # Calculate the distance in the Galaxy's frame
        calculated_distance = v * dt_galaxy
    except ValueError as e:
        return f"Calculation error: {e}"

    # --- Verify the answer ---
    # Check if the proposed answer letter is a valid option
    if proposed_answer_letter not in options:
        return f"Invalid answer: The letter '{proposed_answer_letter}' is not among the options A, B, C, D."

    proposed_answer_value = options[proposed_answer_letter]

    # Find which option is numerically closest to the calculated result
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_distance))

    # Check if the proposed answer is the closest option
    if proposed_answer_letter == closest_option:
        # The answer is correct because it's the best fit among the choices.
        return "Correct"
    else:
        # The answer is incorrect. Provide the reason.
        return (f"Incorrect. The calculated distance is approximately {calculated_distance:.2f} km. "
                f"The proposed answer '{proposed_answer_letter}' ({proposed_answer_value} km) is not the closest option. "
                f"The closest option is '{closest_option}' ({options[closest_option]} km).")

# Run the check
result = check_correctness()
print(result)