import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1 & 2: Define sound sources and the new listener position
    my_new_position = (0, 25)
    
    sources = {
        'dog':          {'pos': (-25, 0), 'Lp_at_1m': 55},
        'train':        {'pos': (50, 0),  'Lp_at_1m': 110},
        'construction': {'pos': (0, 75),  'Lp_at_1m': 90},
        'people':       {'pos': (0, -10), 'Lp_at_1m': 75}
    }

    sum_of_intensities_ratio = 0
    calculated_levels = {}

    # Step 3 & 4: Calculate distance and new sound level for each source
    for name, properties in sources.items():
        source_pos = properties['pos']
        Lp_at_1m = properties['Lp_at_1m']

        # Calculate the new distance from the listener
        distance = math.sqrt((source_pos[0] - my_new_position[0])**2 + (source_pos[1] - my_new_position[1])**2)

        # Calculate the sound level at the new distance
        # Lp = Lp_at_1m - 20 * log10(r)
        Lp_at_new_pos = Lp_at_1m - 20 * math.log10(distance)
        calculated_levels[name] = Lp_at_new_pos
        
        # Add the source's contribution to the total intensity sum
        sum_of_intensities_ratio += 10**(Lp_at_new_pos / 10)

    # Step 5: Convert the total intensity back to decibels
    total_Lp = 10 * math.log10(sum_of_intensities_ratio)

    # Output the final equation and the result
    L_dog = calculated_levels['dog']
    L_train = calculated_levels['train']
    L_construction = calculated_levels['construction']
    L_people = calculated_levels['people']

    print("The sound level (Lp) from each source at the new location is:")
    print(f"  - Dog: {L_dog:.2f} dB")
    print(f"  - Train: {L_train:.2f} dB")
    print(f"  - Construction: {L_construction:.2f} dB")
    print(f"  - People: {L_people:.2f} dB\n")
    
    print("The total sound level is calculated by summing the intensities:")
    print(f"L_total = 10 * log10(10^({L_dog:.2f}/10) + 10^({L_train:.2f}/10) + 10^({L_construction:.2f}/10) + 10^({L_people:.2f}/10))")
    print(f"\nThe total sound level at your new location is: {total_Lp:.2f} dB")
    
    return total_Lp

# Execute the function and store the final answer
final_answer = calculate_total_sound_level()