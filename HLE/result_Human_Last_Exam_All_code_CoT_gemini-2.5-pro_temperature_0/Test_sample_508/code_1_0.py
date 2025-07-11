import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a new location.
    """
    # 1. Initial positions and sound levels at 1 meter
    sources = {
        'dog': {'pos': (-25, 0), 'lp1': 55},
        'train': {'pos': (50, 0), 'lp1': 110},
        'construction': {'pos': (0, 75), 'lp1': 90},
        'people': {'pos': (0, -10), 'lp1': 75}
    }
    
    # 2. Observer's new position
    observer_new_pos = (0, 25)
    
    # Dictionary to store the calculated sound level at the new location for each source
    new_sound_levels = {}
    
    # 3 & 4. Calculate new distances and new sound levels
    for name, data in sources.items():
        source_pos = data['pos']
        lp1 = data['lp1']
        
        # Calculate new distance (r2)
        distance = math.sqrt((source_pos[0] - observer_new_pos[0])**2 + (source_pos[1] - observer_new_pos[1])**2)
        
        # Calculate new sound level (Lp2) at the new distance
        # Lp2 = Lp1 - 20 * log10(r2 / r1), where r1 = 1m
        lp2 = lp1 - 20 * math.log10(distance)
        new_sound_levels[name] = lp2

    # 5. Combine the sound levels
    sum_of_intensities = 0
    for lp in new_sound_levels.values():
        sum_of_intensities += 10**(lp / 10)
        
    total_lp = 10 * math.log10(sum_of_intensities)

    # Print the final equation and result
    lp_dog = new_sound_levels['dog']
    lp_train = new_sound_levels['train']
    lp_construction = new_sound_levels['construction']
    lp_people = new_sound_levels['people']

    equation_str = (
        f"Total Sound Level (dB) = 10 * log10(10^({lp_dog:.2f}/10) + 10^({lp_train:.2f}/10) + "
        f"10^({lp_construction:.2f}/10) + 10^({lp_people:.2f}/10))"
    )
    
    print(equation_str)
    print(f"Total Sound Level (dB) = {total_lp:.2f}")

calculate_total_sound_level()
<<<75.12>>>