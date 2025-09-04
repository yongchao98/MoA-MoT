import math

def check_correctness():
    """
    Checks the correctness of the answer to the special relativity problem.
    """
    # --- Define constants and given values from the question ---
    # Velocity of the ejecta (v) in km/s
    v = 60000
    # Proper time in the ejecta's frame (dt0) in seconds
    dt0 = 50
    # Speed of light (c) in km/s, using the approximation common in these problems
    c = 300000

    # --- The options as provided in the question ---
    # A) 3,060,000 km.
    # B) 2,880,000 km.
    # C) 3,000,000 km.
    # D) 2,940,000 km.
    options = {
        'A': 3060000,
        'B': 2880000,
        'C': 3000000,
        'D': 2940000
    }

    # The final answer to be checked
    final_answer_key = 'A'

    # --- Perform the physics calculation ---
    try:
        # 1. Calculate the ratio v/c and its square
        beta_squared = (v / c)**2

        # 2. Calculate the Lorentz factor (gamma)
        # gamma = 1 / sqrt(1 - v^2/c^2)
        gamma = 1 / math.sqrt(1 - beta_squared)

        # 3. Calculate the time elapsed in the Galaxy's frame (dt) using time dilation
        # dt = gamma * dt0
        dt = gamma * dt0

        # 4. Calculate the distance traveled in the Galaxy's frame
        # distance = v * dt
        calculated_distance = v * dt
    except ValueError as e:
        return f"Calculation error: {e}. This might happen if v >= c."

    # --- Verify the answer ---
    # Find which option is numerically closest to our calculated distance
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_distance))

    # Check if the provided final answer key matches the key of the closest option
    if final_answer_key == closest_option_key:
        # The answer is correct as it identifies the closest option.
        return "Correct"
    else:
        # The answer is incorrect. Provide a detailed reason.
        reason = (
            f"Incorrect. The provided answer is '{final_answer_key}', but the correct answer should be '{closest_option_key}'.\n"
            f"Here is the step-by-step calculation:\n"
            f"1. The velocity (v) is {v} km/s, and the speed of light (c) is {c} km/s. The ratio v/c is {v/c}.\n"
            f"2. The Lorentz factor (gamma) is 1 / sqrt(1 - (v/c)^2) = 1 / sqrt(1 - {beta_squared}) ≈ {gamma:.5f}.\n"
            f"3. The time in the Galaxy's frame (dilated time) is gamma * proper_time = {gamma:.5f} * {dt0} s ≈ {dt:.3f} s.\n"
            f"4. The distance in the Galaxy's frame is v * dilated_time = {v} km/s * {dt:.3f} s ≈ {calculated_distance:.2f} km.\n"
            f"5. The calculated distance {calculated_distance:.2f} km is closest to option {closest_option_key} ({options[closest_option_key]} km), not option {final_answer_key} ({options[final_answer_key]} km)."
        )
        return reason

# Run the check
result = check_correctness()
print(result)