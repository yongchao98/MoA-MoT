import math

def check_supernova_distance():
    """
    This function checks the correctness of the answer to the supernova ejecta problem.

    The problem involves special relativity and time dilation.
    - v: velocity of the ejecta = 60,000 km/s
    - dt0: proper time in the ejecta's frame = 50 s
    - c: speed of light = 300,000 km/s

    The distance in the Galaxy's frame (d) is calculated as d = v * dt,
    where dt is the dilated time in the Galaxy's frame.
    The dilated time is dt = gamma * dt0.
    The Lorentz factor (gamma) is 1 / sqrt(1 - (v/c)^2).

    Therefore, d = v * (1 / sqrt(1 - (v/c)^2)) * dt0.
    """
    # --- Given values and constraints from the question ---
    v = 60000  # km/s
    dt0 = 50   # s (proper time)
    c = 300000 # km/s

    # --- Options from the question ---
    options = {
        'A': 2880000,
        'B': 3000000,
        'C': 2940000,
        'D': 3060000
    }

    # --- The answer to be checked ---
    provided_answer_letter = 'D'
    
    # --- Perform the calculation ---
    try:
        # Calculate the Lorentz factor (gamma)
        beta_squared = (v / c)**2
        if beta_squared >= 1:
            return "Constraint violated: Velocity cannot be greater than or equal to the speed of light."
        gamma = 1 / math.sqrt(1 - beta_squared)

        # Calculate the dilated time in the Galaxy's frame
        dt_galaxy = gamma * dt0

        # Calculate the distance in the Galaxy's frame
        calculated_distance = v * dt_galaxy
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Check the answer ---
    # Find which option is closest to the calculated result
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_distance))

    # Compare the closest option with the provided answer
    if closest_option == provided_answer_letter:
        return "Correct"
    else:
        reason = (
            f"The provided answer is incorrect.\n"
            f"The calculated distance is approximately {calculated_distance:,.2f} km.\n"
            f"The non-relativistic distance (ignoring time dilation) would be {v * dt0:,} km, which is option B.\n"
            f"Due to time dilation, the distance must be greater than 3,000,000 km.\n"
            f"The calculated value {calculated_distance:,.2f} km is closest to option {closest_option} ({options[closest_option]:,} km).\n"
            f"The provided answer was {provided_answer_letter}, but the correct option is {closest_option}."
        )
        return reason

# Run the check and print the result
print(check_supernova_distance())