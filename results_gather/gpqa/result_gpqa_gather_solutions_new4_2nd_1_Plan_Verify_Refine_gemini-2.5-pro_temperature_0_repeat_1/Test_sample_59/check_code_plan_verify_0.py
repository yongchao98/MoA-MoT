import math

def check_correctness():
    """
    Checks the correctness of the answer to the special relativity problem.

    The problem involves calculating the proper time (time experienced by the astronaut)
    for a journey from the Large Magellanic Cloud (LMC) to Earth.

    The solution requires:
    1. Using an external value for the distance to the LMC.
    2. Calculating the Lorentz factor (gamma).
    3. Calculating the time in Earth's frame of reference.
    4. Calculating the proper time for the astronaut using the time dilation formula.
    5. Checking if the astronaut survives the journey.
    """

    # --- Problem Parameters ---
    v_ratio = 0.99999987  # v/c
    astronaut_start_age = 22
    alien_lifespan = 150

    # --- Answer to Verify ---
    # The provided answer is C, which corresponds to 81 years.
    # This implies the astronaut survives.
    expected_time_years = 81
    expected_survival = True

    # --- External Data ---
    # The distance to the LMC is not given. The analysis correctly assumes a standard
    # value. We will use a range of plausible values to ensure robustness.
    # Common values range from ~158,000 to ~163,000 light-years.
    plausible_distances_ly = [159000, 160000, 163000]
    
    # --- Calculation ---
    try:
        # Calculate the Lorentz factor (gamma), which only depends on velocity.
        gamma = 1 / math.sqrt(1 - v_ratio**2)
    except ValueError:
        return "Calculation Error: The velocity ratio v/c cannot be 1 or greater."

    # Check if any plausible distance yields a result consistent with the answer.
    for distance_ly in plausible_distances_ly:
        # Calculate time as measured from Earth (delta_t)
        delta_t_earth = distance_ly / v_ratio

        # Calculate proper time for the astronaut (delta_t_0)
        delta_t_astronaut = delta_t_earth / gamma

        # Check survival condition
        final_age = astronaut_start_age + delta_t_astronaut
        astronaut_survives = final_age < alien_lifespan

        # Check if this calculation matches the expected answer
        # We use a tolerance because the distance is an approximation and the answer is rounded.
        # A 2% relative tolerance is reasonable.
        time_matches = math.isclose(delta_t_astronaut, expected_time_years, rel_tol=0.02)
        
        if astronaut_survives == expected_survival and time_matches:
            return "Correct"

    # If the loop finishes without finding a match, the answer is incorrect.
    # We perform one last calculation with the most common distance (160,000 ly) for the error message.
    distance_ly = 160000
    delta_t_earth = distance_ly / v_ratio
    delta_t_astronaut = delta_t_earth / gamma
    final_age = astronaut_start_age + delta_t_astronaut
    astronaut_survives = final_age < alien_lifespan

    if astronaut_survives != expected_survival:
        return (f"Incorrect. The survival condition is not met. "
                f"The answer implies the astronaut survives, but the calculation shows their final age "
                f"would be {final_age:.2f} years, which is {'not ' if astronaut_survives else ''}less than the {alien_lifespan}-year lifespan.")

    return (f"Incorrect. The calculated travel time does not match the answer. "
            f"Using a standard distance of {distance_ly} light-years, the calculated travel time for the astronaut is {delta_t_astronaut:.2f} years. "
            f"This is not sufficiently close to the given answer of {expected_time_years} years. "
            f"Even checking other plausible distances for the LMC (from 159k to 163k light-years) yields results clustered around 81-83 years, but none match the options other than {expected_time_years} which is confirmed by the calculation.")

# The code above is a function definition. To get the output, we need to call it.
# In a real execution environment, this would be run.
# For this format, we will print the result of the function call.
# print(check_correctness())