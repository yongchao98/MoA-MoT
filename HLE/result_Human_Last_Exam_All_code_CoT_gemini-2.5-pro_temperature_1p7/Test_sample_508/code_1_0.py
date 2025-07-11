import math

def calculate_total_sound_level():
    """
    This script calculates the total sound level experienced by an observer
    from multiple sound sources at various distances.
    """
    # 1. Define sound sources with their sound pressure level (SPL) at 1m
    #    and their coordinates.
    sources = [
        {'name': 'Dog', 'spl_1m': 55, 'pos': (-25, 0)},
        {'name': 'Train', 'spl_1m': 110, 'pos': (50, 0)},
        {'name': 'Construction', 'spl_1m': 90, 'pos': (0, 75)},
        {'name': 'People', 'spl_1m': 75, 'pos': (0, -10)}
    ]

    # 2. Define the observer's new position.
    # The observer starts at (0,0) and walks 25m towards the construction at (0, 75).
    my_new_pos = (0, 25)
    reference_distance = 1.0  # All source SPLs are measured at 1 meter.

    print("--- Individual Sound Level Calculations ---")

    individual_levels = {}
    sum_of_intensities = 0

    # Loop through each source to calculate its contribution
    for source in sources:
        # 3. Calculate the new distance from the observer to the source.
        distance_to_source = math.sqrt(
            (source['pos'][0] - my_new_pos[0])**2 + 
            (source['pos'][1] - my_new_pos[1])**2
        )

        # 4. Calculate the attenuated sound level at the observer's location.
        # Formula: Lp2 = Lp1 - 20 * log10(r2 / r1)
        spl_at_location = source['spl_1m'] - 20 * math.log10(distance_to_source / reference_distance)
        individual_levels[source['name']] = spl_at_location
        
        # 5. Convert the dB level to a relative intensity and add it to a running sum.
        # Formula: intensity = 10^(Lp/10)
        sum_of_intensities += 10**(spl_at_location / 10)

        print(f"Calculated level from {source['name']}: {spl_at_location:.2f} dB at {distance_to_source:.2f} meters")

    # 6. Combine the sound levels by converting the total intensity back to dB.
    # Formula: L_total = 10 * log10(sum of intensities)
    total_spl = 10 * math.log10(sum_of_intensities)

    print("\n--- Total Sound Level Calculation ---")

    # Retrieve individual calculated SPLs for the final equation display
    dog_spl = individual_levels['Dog']
    train_spl = individual_levels['Train']
    construction_spl = individual_levels['Construction']
    people_spl = individual_levels['People']

    # Display the final equation with all numbers plugged in, as requested.
    equation_str = (
        f"Total dB = 10 * log10(10^({dog_spl:.2f}/10) + 10^({train_spl:.2f}/10) + "
        f"10^({construction_spl:.2f}/10) + 10^({people_spl:.2f}/10))"
    )

    print("Final Equation:")
    print(equation_str)

    print(f"\nThe total sound level you hear at your new location is: {total_spl:.2f} dB")
    
    # Return the final numeric answer for wrapping.
    return total_spl

# Execute the function and capture the result
final_answer = calculate_total_sound_level()