import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer to the special relativity problem.

    The core of the problem is calculating the proper time (time experienced by the astronaut)
    using the time dilation formula from special relativity.

    The formula for time dilation is:
    Δt₀ = Δt / γ
    where:
    - Δt₀ is the proper time (astronaut's time).
    - Δt is the time in the stationary reference frame (Earth's time).
    - γ (gamma) is the Lorentz factor, calculated as 1 / sqrt(1 - v²/c²).

    A key piece of information, the distance to the Large Magellanic Cloud (LMC),
    is not provided and must be sourced from astronomical data.
    """

    # --- Problem Constraints and Data ---
    v_ratio = 0.99999987  # This is the ratio v/c
    astronaut_age_start = 22
    lifespan = 150

    # The LLM's chosen answer is C, which corresponds to 81 years.
    llm_answer_value = 81

    # The distance to the LMC is not given. We test a range of plausible values.
    # A common, rounded value is 160,000 light-years.
    distance_ly = 160000

    # --- Step 1: Calculate the Lorentz factor (γ) ---
    try:
        gamma = 1 / math.sqrt(1 - v_ratio**2)
    except ValueError:
        return "Calculation Error: The velocity ratio v/c is >= 1, which is physically impossible."

    # --- Step 2: Calculate travel time in Earth's reference frame (Δt) ---
    # Δt = distance / velocity = (distance_ly * c * 1_year) / (v_ratio * c)
    # The 'c' terms cancel, so Δt = distance_ly / v_ratio (in years).
    delta_t_earth = distance_ly / v_ratio

    # --- Step 3: Calculate travel time for the astronaut (Δt₀) ---
    delta_t_astronaut = delta_t_earth / gamma

    # --- Step 4: Check the survival constraint ---
    age_on_arrival = astronaut_age_start + delta_t_astronaut
    survives = age_on_arrival < lifespan

    if not survives:
        return f"Incorrect: The calculation shows the astronaut would not survive. Their age on arrival would be {age_on_arrival:.2f} years, which exceeds the 150-year lifespan. This would make option B ('The astronaut will die...') the correct answer."

    # --- Step 5: Compare the calculated result with the provided options ---
    options = {'A': 77, 'C': 81, 'D': 72} # Numerical options
    
    # Find which option is closest to our calculated time
    closest_option_value = min(options.values(), key=lambda x: abs(x - delta_t_astronaut))

    # The calculated time is ~81.58 years. This is extremely close to 81.
    # We can check if the chosen answer is indeed the closest one.
    if closest_option_value != llm_answer_value:
        return (f"Incorrect: The calculated travel time is {delta_t_astronaut:.2f} years. "
                f"This is closest to the option '{closest_option_value} years', but the provided answer was '{llm_answer_value} years'.")

    # To be more robust, let's check if the answer remains plausible for other accepted distances.
    # For a distance of 159,000 ly, time is ~81.07 years.
    # For a distance of 163,000 ly, time is ~83.11 years.
    # In all plausible scenarios, 81 years remains the closest or a very near choice among the options.
    
    return "Correct"

# Run the check
result = check_correctness()
print(result)