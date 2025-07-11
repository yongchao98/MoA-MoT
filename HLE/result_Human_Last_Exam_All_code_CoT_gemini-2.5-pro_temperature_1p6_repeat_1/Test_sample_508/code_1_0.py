import math

def calculate_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial sound source data
    sources = {
        'dog': {'L1': 55, 'pos': (-25, 0)},
        'train': {'L1': 110, 'pos': (50, 0)},
        'construction': {'L1': 90, 'pos': (0, 75)},
        'people': {'L1': 75, 'pos': (0, -10)}
    }

    # 2. Define the observer's new position
    my_pos = (0, 25)

    sound_levels_at_new_pos = {}
    sum_of_intensities = 0

    print("Calculating individual sound levels at the new location:\n")

    for name, data in sources.items():
        # 3. Calculate new distance to the source
        source_pos = data['pos']
        distance = math.sqrt((source_pos[0] - my_pos[0])**2 + (source_pos[1] - my_pos[1])**2)

        # 4. Calculate the sound level (L2) from this source at the new distance
        l1 = data['L1']
        # The formula is L2 = L1 - 20 * log10(r2 / r1). Here r1=1, so it simplifies.
        l2 = l1 - 20 * math.log10(distance)
        
        sound_levels_at_new_pos[name] = l2
        
        print(f"Sound level from {name}:")
        print(f"  - Distance from new position: {distance:.2f} meters")
        print(f"  - Calculated sound level: {l2:.2f} dB")
        
        # Prepare for total calculation
        sum_of_intensities += 10**(l2 / 10)

    # 5. Calculate the total sound level
    total_sound_level = 10 * math.log10(sum_of_intensities)

    L_dog = sound_levels_at_new_pos['dog']
    L_train = sound_levels_at_new_pos['train']
    L_construction = sound_levels_at_new_pos['construction']
    L_people = sound_levels_at_new_pos['people']
    
    print("\n----------------------------------------------------")
    print("Calculating the total sound level:")
    print("Formula: L_total = 10 * log10(10^(L1/10) + 10^(L2/10) + ...)")
    print(f"L_total = 10 * log10(10^({L_dog:.2f}/10) + 10^({L_train:.2f}/10) + 10^({L_construction:.2f}/10) + 10^({L_people:.2f}/10))")
    print(f"The total sound level you hear is: {total_sound_level:.2f} dB")
    print("----------------------------------------------------\n")

    return total_sound_level

if __name__ == '__main__':
    final_answer = calculate_sound_level()
    print(f'<<<{final_answer:.2f}>>>')