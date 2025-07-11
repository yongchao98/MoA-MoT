import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1: Define initial data for each sound source
    sources = {
        'Dog':          {'pos': (-25, 0), 'L_1m': 55},
        'Train':        {'pos': (50, 0),  'L_1m': 110},
        'Construction': {'pos': (0, 75),  'L_1m': 90},
        'People':       {'pos': (0, -10), 'L_1m': 75}
    }

    # Step 2: Define your new position
    my_new_pos = (0, 25)

    new_levels_db = {}
    intensity_sum_str_parts = []
    total_intensity = 0

    print("Step 1: Calculate the sound level (dB) from each source at your new location.\n")
    # Step 3 & 4: Calculate for each source
    for name, data in sources.items():
        source_pos = data['pos']
        L_1m = data['L_1m']

        # Calculate the distance from your new position to the sound source
        distance = math.sqrt((source_pos[0] - my_new_pos[0])**2 + (source_pos[1] - my_new_pos[1])**2)

        # Calculate the new sound level (L_new) at your location using the attenuation formula
        L_new = L_1m - 20 * math.log10(distance)
        new_levels_db[name] = L_new

        # Convert the new dB level to a relative intensity and add it to the total
        intensity = 10**(L_new / 10)
        total_intensity += intensity

        # Prepare strings for the final equation output
        intensity_sum_str_parts.append(f"10^({L_new:.2f}/10)")
        
        print(f"Source: {name}")
        print(f"  - Distance to new location: {distance:.2f} m")
        print(f"  - Sound level at new location: {L_new:.2f} dB\n")

    # Convert the total intensity back to the final decibel level
    L_total = 10 * math.log10(total_intensity)

    # Construct the final equation string showing all the numbers
    equation_str = " + ".join(intensity_sum_str_parts)

    print("Step 2: Combine the sound levels to find the total.")
    print("The total sound level is calculated by summing the intensities and converting back to decibels.\n")
    print("Final Equation:")
    print(f"Total Sound Level = 10 * log10({equation_str})")
    print(f"\nResult:")
    print(f"Total Sound Level = {L_total:.2f} dB")


if __name__ == '__main__':
    calculate_total_sound_level()
