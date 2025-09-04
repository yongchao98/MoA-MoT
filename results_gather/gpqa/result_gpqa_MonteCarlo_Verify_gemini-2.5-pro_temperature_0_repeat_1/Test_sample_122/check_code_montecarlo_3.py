import math

def check_answer():
    """
    Checks the correctness of the given answer for the special relativity problem.
    """
    # --- Problem Constants and Given Values ---
    # Speed of light in km/s
    c = 299792.458
    # Relative velocity of the ejecta in km/s
    v = 60000
    # Proper time elapsed in the ejecta's frame in seconds
    dt_proper = 50

    # --- The LLM's Answer ---
    # The chosen option is B
    llm_choice_key = "B"
    options = {
        "A": 2940000,
        "B": 3060000,
        "C": 3000000,
        "D": 2880000
    }
    llm_choice_value = options[llm_choice_key]

    # --- Step 1: Verify the logical estimation/sanity check ---
    # Non-relativistic distance calculation
    d_non_relativistic = v * dt_proper
    if not math.isclose(d_non_relativistic, 3000000):
        return f"Constraint check failed: The non-relativistic distance should be v * t = 60000 * 50 = 3,000,000 km. The LLM's reasoning uses this correct value."

    # The core of the logical estimation is that due to time dilation, the distance
    # must be greater than the non-relativistic distance.
    if not (llm_choice_value > d_non_relativistic):
        return f"Constraint check failed: The answer must be greater than the non-relativistic distance of {d_non_relativistic} km. The chosen answer {llm_choice_value} km does not satisfy this."

    # Check if B is the only option greater than 3,000,000 km
    plausible_options = [k for k, val in options.items() if val > d_non_relativistic]
    if plausible_options != ['B']:
        return f"Constraint check failed: The LLM's reasoning that only option B is greater than 3,000,000 km is incorrect. The plausible options are {plausible_options}."

    # --- Step 2: Perform the precise calculation ---
    try:
        # Calculate the Lorentz factor (gamma)
        beta_squared = (v / c)**2
        gamma = 1 / math.sqrt(1 - beta_squared)

        # Calculate the time elapsed in the Galaxy frame (time dilation)
        dt_galaxy = gamma * dt_proper

        # Calculate the distance traveled in the Galaxy frame
        d_galaxy_calculated = v * dt_galaxy
    except Exception as e:
        return f"An error occurred during the physics calculation: {e}"

    # --- Step 3: Verify the final answer ---
    # Find the option that is numerically closest to our precise calculation.
    closest_option_key = min(options.keys(), key=lambda k: abs(options[k] - d_galaxy_calculated))

    # Check if the LLM's chosen option is the closest one.
    if llm_choice_key != closest_option_key:
        return (f"Incorrect. The calculated distance is {d_galaxy_calculated:.2f} km. "
                f"The closest option is {closest_option_key} ({options[closest_option_key]} km), "
                f"but the provided answer chose {llm_choice_key} ({llm_choice_value} km).")

    # If all checks pass, the reasoning and the final choice are correct.
    return "Correct"

# Execute the check
result = check_answer()
if result == "Correct":
    print("Correct")
else:
    print(result)