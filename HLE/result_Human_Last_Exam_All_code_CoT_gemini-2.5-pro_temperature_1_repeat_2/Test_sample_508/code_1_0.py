import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial positions and sound levels (at 1 meter)
    sources = {
        'dog': {'pos': (-25, 0), 'L0': 55},
        'train': {'pos': (50, 0), 'L0': 110},
        'construction': {'pos': (0, 75), 'L0': 90},
        'people': {'pos': (0, -10), 'L0': 75}
    }

    # 2. Define the observer's new position
    observer_pos = (0, 25)

    # Dictionary to store the calculated sound levels at the new position
    new_sound_levels = {}
    
    # 3 & 4. Calculate new distances and the sound level from each source at the new location
    for name, data in sources.items():
        source_pos = data['pos']
        L0 = data['L0']
        
        # Calculate distance
        distance = math.sqrt((source_pos[0] - observer_pos[0])**2 + (source_pos[1] - observer_pos[1])**2)
        
        # Calculate sound level at the new distance. Avoid log10(0) or negative distances.
        if distance > 0:
            L = L0 - 20 * math.log10(distance)
            new_sound_levels[name] = L
        else:
            # This case shouldn't happen with the given problem data
            new_sound_levels[name] = float('inf')

    l_dog = new_sound_levels['dog']
    l_train = new_sound_levels['train']
    l_construction = new_sound_levels['construction']
    l_people = new_sound_levels['people']

    # 5. Calculate the total sound level by summing intensities
    total_intensity = (10**(l_dog / 10) + 
                       10**(l_train / 10) + 
                       10**(l_construction / 10) + 
                       10**(l_people / 10))
    
    total_sound_level = 10 * math.log10(total_intensity)

    # Print the results and the final equation
    print(f"After moving, the individual sound levels are:")
    print(f"Dog: {l_dog:.2f} dB")
    print(f"Train: {l_train:.2f} dB")
    print(f"Construction: {l_construction:.2f} dB")
    print(f"People: {l_people:.2f} dB\n")

    print("The total sound level is calculated as follows:")
    equation = (
        f"Total dB = 10 * log10(10^({l_dog:.2f}/10) + "
        f"10^({l_train:.2f}/10) + "
        f"10^({l_construction:.2f}/10) + "
        f"10^({l_people:.2f}/10))"
    )
    print(equation)
    print(f"\nTotal Sound Level = {total_sound_level:.2f} dB")
    
    return total_sound_level

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    print(f"\n<<<{final_answer:.2f}>>>")