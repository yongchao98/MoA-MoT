import math

def calculate_total_sound_level():
    """
    This function calculates the total sound level from multiple sources at a new location.
    """
    # Step 1 & 2: Define sound sources and the observer's new position.
    # Your new position is (0, 25) after walking 25m towards the construction site from the origin.
    my_pos = (0, 25)

    sources = [
        {'name': 'dog', 'level_at_1m': 55, 'pos': (-25, 0)},
        {'name': 'train', 'level_at_1m': 110, 'pos': (50, 0)},
        {'name': 'construction', 'level_at_1m': 90, 'pos': (0, 75)},
        {'name': 'people', 'level_at_1m': 75, 'pos': (0, -10)}
    ]

    new_levels = []
    intensity_terms = []
    
    print("Calculating the sound level from each source at the new position (0, 25)...")

    # Step 3 & 4: Calculate new distance and new sound level for each source.
    for source in sources:
        distance = math.sqrt((source['pos'][0] - my_pos[0])**2 + (source['pos'][1] - my_pos[1])**2)
        # Formula: L_new = L_at_1m - 20 * log10(distance)
        level_at_new_pos = source['level_at_1m'] - 20 * math.log10(distance)
        new_levels.append(level_at_new_pos)
        print(f"  - {source['name'].capitalize()}: Level at {distance:.2f}m is {level_at_new_pos:.2f} dB")

    # Step 5: Convert levels to intensities, sum them, and convert back to total dB.
    total_intensity = sum(10**(level / 10) for level in new_levels)
    total_sound_level = 10 * math.log10(total_intensity)

    print("\nTo find the total sound level, we combine the sources using their intensities:")
    
    # Building and printing the equation with all the numbers
    equation_levels_str = " + ".join([f"10^({level:.2f}/10)" for level in new_levels])
    print(f"L_total = 10 * log10({equation_levels_str})")
    
    intensities = [10**(level / 10) for level in new_levels]
    equation_intensities_str = " + ".join([f"{intensity:.2f}" for intensity in intensities])
    print(f"L_total = 10 * log10({equation_intensities_str})")

    print(f"L_total = 10 * log10({total_intensity:.2f})")
    print(f"L_total = {total_sound_level:.2f} dB")
    
    print("\n-------------------------------------------------")
    print(f"The total sound level I hear is {total_sound_level:.2f} dB.")
    print("-------------------------------------------------")
    
    return total_sound_level

if __name__ == "__main__":
    final_answer = calculate_total_sound_level()
    print(f"<<<{final_answer:.2f}>>>")
