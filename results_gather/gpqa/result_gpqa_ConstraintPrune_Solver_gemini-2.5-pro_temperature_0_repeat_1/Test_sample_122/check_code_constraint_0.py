import math

def check_supernova_distance():
    """
    This function verifies the answer to the special relativity problem.
    It calculates the distance traveled by the ejecta in the Galaxy's reference frame
    and compares it with the provided answer.
    """

    # --- Problem Parameters and Constants ---
    # Velocity of the ejecta relative to the Galaxy (v)
    v = 60000  # km/s
    # Time elapsed in the ejecta's reference frame (proper time, Δt')
    dt_prime = 50  # s
    # Speed of light in a vacuum (c)
    c = 299792.458  # km/s

    # --- Provided Answer ---
    # The LLM's answer is B, which corresponds to 3,060,000 km.
    provided_answer_label = "B"
    options = {
        "A": 2880000,
        "B": 3060000,
        "C": 3000000,
        "D": 2940000
    }
    provided_answer_value = options[provided_answer_label]

    # --- Step 1: Check the non-relativistic case (a common error) ---
    # This is the distance if we ignore time dilation.
    # d_non_relativistic = v * Δt'
    d_non_relativistic = v * dt_prime
    
    # The answer correctly states that the relativistic distance must be greater than this.
    if provided_answer_value <= d_non_relativistic:
        return (f"Incorrect. The provided answer {provided_answer_value:,} km is less than or equal to "
                f"the non-relativistic distance of {d_non_relativistic:,} km. Due to time dilation, "
                "the time elapsed in the Galaxy frame is longer, so the distance traveled must be greater.")

    # --- Step 2: Calculate the correct relativistic distance ---
    # Calculate the Lorentz factor (gamma)
    # gamma = 1 / sqrt(1 - v^2/c^2)
    try:
        beta_squared = (v / c)**2
        gamma = 1 / math.sqrt(1 - beta_squared)
    except ValueError:
        return "Incorrect. Calculation failed because velocity is greater than or equal to the speed of light."

    # Calculate the time elapsed in the Galaxy's reference frame (Δt)
    # Δt = gamma * Δt'
    dt_galaxy = gamma * dt_prime

    # Calculate the distance traveled in the Galaxy's reference frame (d)
    # d = v * Δt
    d_galaxy_calculated = v * dt_galaxy

    # --- Step 3: Compare the calculated distance with the provided answer ---
    # For multiple-choice questions, the answer might be rounded.
    # We check if the provided answer is the closest option to our calculated value.

    # Find the option closest to our calculated result
    closest_option_label = min(options, key=lambda k: abs(options[k] - d_galaxy_calculated))

    if closest_option_label != provided_answer_label:
        return (f"Incorrect. The calculated distance is approximately {d_galaxy_calculated:,.2f} km. "
                f"The closest option is {closest_option_label} ({options[closest_option_label]:,} km), "
                f"not {provided_answer_label} ({provided_answer_value:,} km).")

    # Check if the value of the closest option is reasonably close to the calculated value.
    # A tolerance of 1% is reasonable for such problems.
    tolerance = 0.01
    difference_ratio = abs(provided_answer_value - d_galaxy_calculated) / d_galaxy_calculated

    if difference_ratio > tolerance:
        return (f"Incorrect. Although option {provided_answer_label} is the closest, its value ({provided_answer_value:,} km) "
                f"differs from the calculated value ({d_galaxy_calculated:,.2f} km) by more than {tolerance:.1%}. "
                f"The discrepancy is {difference_ratio:.2%}, which may indicate an error in the question's options or rounding.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_supernova_distance()
print(result)