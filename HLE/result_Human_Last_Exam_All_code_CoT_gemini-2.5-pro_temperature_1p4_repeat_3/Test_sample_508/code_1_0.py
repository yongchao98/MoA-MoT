import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # My new position after walking 25m towards the construction
    my_pos = (0, 25)

    # Sound sources: {name: (level_at_1m, (x, y))}
    sources = {
        'dog':          (55, (-25, 0)),
        'train':        (110, (50, 0)),
        'construction': (90, (0, 75)),
        'people':       (75, (0, -10)),
    }

    total_relative_intensity = 0
    new_sound_levels = {}
    
    print("Calculating sound level from each source at your new location...\n")

    # Calculate sound level and relative intensity for each source
    for name, (level_at_1m, pos) in sources.items():
        # Calculate distance from my new position to the source
        distance = math.sqrt((my_pos[0] - pos[0])**2 + (my_pos[1] - pos[1])**2)
        
        # Calculate the sound level at the new distance (attenuation)
        level_at_distance = level_at_1m - 20 * math.log10(distance)
        new_sound_levels[name] = level_at_distance
        
        # Convert sound level (dB) to relative intensity
        relative_intensity = 10**(level_at_distance / 10)
        
        # Add to the total intensity
        total_relative_intensity += relative_intensity
        
        print(f"The distance to the {name} is {distance:.2f} meters.")
        print(f"The sound level from the {name} is {level_at_distance:.2f} dB.\n")

    # Convert the total relative intensity back to decibels
    total_sound_level = 10 * math.log10(total_relative_intensity)

    # Print the final combined calculation
    # Building the string for the equation
    equation_parts = []
    for level in new_sound_levels.values():
        equation_parts.append(f"10^({level:.2f}/10)")
    equation_str = " + ".join(equation_parts)
    
    print("--------------------------------------------------")
    print("The total sound level is calculated by converting each dB value to intensity,")
    print("summing the intensities, and converting back to dB.\n")
    print("Final Equation:")
    print(f"Total dB = 10 * log10({equation_str})")
    print(f"Total dB = {total_sound_level:.2f}")

    return total_sound_level

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    # The final answer is wrapped for direct parsing.
    print(f"\n<<<The total sound level you hear is approximately {final_answer:.2f} dB.>>>")
