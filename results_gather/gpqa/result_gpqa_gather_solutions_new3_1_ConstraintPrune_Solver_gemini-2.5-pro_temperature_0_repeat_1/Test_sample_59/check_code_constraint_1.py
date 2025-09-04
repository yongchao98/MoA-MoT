import math

def check_answer():
    """
    Checks the correctness of the answer to the special relativity problem.
    """
    # --- Problem Constraints and Given Data ---
    speed_ratio = 0.99999987  # v/c
    astronaut_initial_age = 22
    alien_lifespan = 150
    
    # The distance to the Large Magellanic Cloud (LMC) is not given,
    # but astronomical estimates range from ~159,000 to ~163,000 light-years.
    # We will test a few values in this range.
    distances_ly = {
        "low_estimate": 159000,
        "common_estimate": 160000,
        "high_estimate": 163000
    }
    
    # --- Options from the Question ---
    # A) 72 years, B) 77 years, C) 81 years, D) The astronaut will die
    options = {'A': 72, 'B': 77, 'C': 81}
    final_answer_choice = 'C'

    # --- Physics Calculations ---
    try:
        # Calculate the Lorentz factor (gamma)
        # gamma = 1 / sqrt(1 - v^2/c^2)
        gamma = 1 / math.sqrt(1 - speed_ratio**2)

        results = {}
        for name, dist in distances_ly.items():
            # Time in Earth's reference frame (t_earth)
            # t_earth = distance / speed
            t_earth = dist / speed_ratio
            
            # Time for the astronaut (proper time, t_astronaut)
            # t_astronaut = t_earth / gamma
            t_astronaut = t_earth / gamma
            results[name] = t_astronaut

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    
    # 1. Check the survival constraint
    # We use the highest calculated travel time for the most stringent check.
    max_travel_time = max(results.values())
    age_on_arrival = astronaut_initial_age + max_travel_time
    
    if age_on_arrival >= alien_lifespan:
        # This would mean option D is correct.
        if final_answer_choice == 'D':
            return "Correct"
        else:
            return f"Incorrect. The calculation shows the astronaut would not survive (age on arrival: {age_on_arrival:.2f} years), but the answer was not D."

    # If the astronaut survives, D is incorrect.
    if final_answer_choice == 'D':
        return f"Incorrect. The calculation shows the astronaut survives the journey (age on arrival: {age_on_arrival:.2f} years), but the answer was D (will die)."

    # 2. Check which numerical option is closest
    # All calculated values should point to the same option.
    is_consistent = True
    target_option = None

    for name, t_astro in results.items():
        # Find the closest option for this calculated time
        closest_option = min(options.keys(), key=lambda k: abs(options[k] - t_astro))
        
        if target_option is None:
            target_option = closest_option
        elif target_option != closest_option:
            is_consistent = False
            break
            
    if not is_consistent:
        return "Incorrect. The choice of distance to LMC leads to different best-fit answers, making the question ambiguous."

    # 3. Final check against the provided answer
    if target_option == final_answer_choice:
        return "Correct"
    else:
        calculated_time = results["common_estimate"]
        return (f"Incorrect. The provided answer is {final_answer_choice} ({options.get(final_answer_choice)} years). "
                f"However, the calculation consistently points to option {target_option} ({options.get(target_option)} years). "
                f"For a distance of 160,000 light-years, the travel time is ~{calculated_time:.2f} years, which is closest to {options.get(target_option)}.")

# Run the check
print(check_answer())