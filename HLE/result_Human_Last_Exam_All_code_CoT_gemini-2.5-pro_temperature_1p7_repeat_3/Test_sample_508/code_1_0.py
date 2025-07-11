import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """

    # 1. Define initial sources and your new position
    sources = {
        'dog': {'L1': 55, 'pos': (-25, 0)},
        'train': {'L1': 110, 'pos': (50, 0)},
        'construction': {'L1': 90, 'pos': (0, 75)},
        'people': {'L1': 75, 'pos': (0, -10)}
    }
    my_new_pos = (0, 25)

    # Dictionary to store the calculated sound level for each source at the new location
    new_sound_levels = {}
    
    # 2. Calculate new sound levels for each source
    for name, data in sources.items():
        # Calculate the new distance from the source to your new position
        distance = math.sqrt((data['pos'][0] - my_new_pos[0])**2 + (data['pos'][1] - my_new_pos[1])**2)
        
        # Calculate the sound level (L2) at the new distance. L1 is the level at 1m.
        # Formula: L2 = L1 - 20 * log10(d2 / d1), where d1 = 1 meter
        if distance > 0:
            l2 = data['L1'] - 20 * math.log10(distance)
            new_sound_levels[name] = l2
        else:
            # Handle the case where you are at the source, which is effectively infinite dB.
            # In a real-world scenario, this is capped. For this problem, it's not possible.
            new_sound_levels[name] = float('inf')


    # 3. Convert dB to intensity, sum them up
    total_intensity = 0
    for level in new_sound_levels.values():
        total_intensity += 10**(level / 10)

    # 4. Convert total intensity back to dB
    # L_total = 10 * log10(Sum of 10^(Li/10))
    if total_intensity > 0:
        total_db = 10 * math.log10(total_intensity)
    else:
        total_db = -float('inf') # Log of 0 is undefined, represents silence

    # 5. Print the final equation and result
    print("The individual sound levels at your new location are:")
    for name, level in new_sound_levels.items():
        print(f"- {name.capitalize()}: {level:.2f} dB")
    
    print("\nThe total sound level is found by summing the sound intensities and converting back to decibels.")
    print("The final calculation is:")
    
    # Construct the equation string
    equation_str = "10 * log10("
    for level in new_sound_levels.values():
        equation_str += f"10^({level:.2f}/10) + "
    # Remove the trailing ' + ' and add the closing parenthesis
    equation_str = equation_str[:-3] + ")"
    
    print(f"{equation_str} = {total_db:.2f} dB")


calculate_total_sound_level()
<<<75.12>>>