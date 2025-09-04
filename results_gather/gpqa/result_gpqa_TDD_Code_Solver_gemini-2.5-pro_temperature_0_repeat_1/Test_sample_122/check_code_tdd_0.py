import math

def check_supernova_answer():
    """
    Checks the correctness of the provided answer for the supernova distance problem.
    """
    # --- Problem Parameters ---
    # Relative velocity in km/s
    v = 60000.0
    # Proper time in the ejecta's frame in seconds
    t_proper = 50.0
    # Speed of light in km/s (standard value)
    c = 299792.458

    # --- Provided Answer ---
    # The LLM's answer is 'A', which corresponds to 3,060,000 km.
    expected_answer_label = 'A'
    options = {
        'A': 3060000.0,
        'B': 2880000.0,
        'C': 3000000.0,
        'D': 2940000.0
    }

    # --- Independent Calculation ---
    # Step 1: Check for non-relativistic conditions
    if abs(v) >= c:
        return "Error: Velocity cannot be equal to or greater than the speed of light."

    # Step 2: Calculate the Lorentz factor (gamma)
    try:
        beta_squared = (v / c) ** 2
        gamma = 1 / math.sqrt(1 - beta_squared)
    except ValueError:
        return "Error during calculation: Invalid value for gamma (velocity might be too high)."

    # Step 3: Calculate the time in the Galaxy's reference frame (time dilation)
    t_galaxy = gamma * t_proper

    # Step 4: Calculate the distance in the Galaxy's reference frame
    calculated_distance = v * t_galaxy

    # --- Verification ---
    # Step 5: Find the closest option to the calculated distance
    closest_option_label = None
    min_difference = float('inf')

    for label, value in options.items():
        difference = abs(calculated_distance - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_label = label

    # Step 6: Check if the closest option matches the provided answer
    if closest_option_label == expected_answer_label:
        # The logic is correct, and the chosen option is the best fit.
        return "Correct"
    else:
        # The logic or the chosen option is incorrect.
        return (f"Incorrect. The calculated distance is approximately {calculated_distance:,.2f} km. "
                f"The closest option is {closest_option_label} ({options[closest_option_label]:,.0f} km), "
                f"but the provided answer was {expected_answer_label} ({options[expected_answer_label]:,.0f} km).")

# Run the check
result = check_supernova_answer()
print(result)