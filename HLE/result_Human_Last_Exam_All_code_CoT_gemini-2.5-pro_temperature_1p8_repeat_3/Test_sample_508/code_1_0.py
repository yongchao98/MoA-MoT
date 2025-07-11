import math

def calculate_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial positions and sound levels (at 1m)
    sources = {
        'dog':          {'pos': (-25, 0), 'L1': 55},
        'train':        {'pos': (50, 0),  'L1': 110},
        'construction': {'pos': (0, 75),  'L1': 90},
        'people':       {'pos': (0, -10), 'L1': 75}
    }

    # 2. Define the listener's new position
    my_new_pos = (0, 25)

    individual_levels = []
    
    # 3. Calculate new distance and new sound level for each source
    for source_name, data in sources.items():
        # Euclidean distance formula: sqrt((x2 - x1)^2 + (y2 - y1)^2)
        distance = math.sqrt((my_new_pos[0] - data['pos'][0])**2 + (my_new_pos[1] - data['pos'][1])**2)
        
        # Sound level formula: L2 = L1 - 20 * log10(r2) (since r1 = 1)
        # L1 is the sound level at 1 meter
        l2 = data['L1'] - 20 * math.log10(distance)
        individual_levels.append(l2)
        # print(f"Sound level from {source_name} at new location: {l2:.2f} dB (distance: {distance:.2f}m)")

    # 4. Combine the sound levels
    # Convert each dB to intensity, sum them, and convert back to dB
    # L_total = 10 * log10(sum(10^(Li/10)))
    total_intensity = sum([10**(level / 10) for level in individual_levels])
    total_level_db = 10 * math.log10(total_intensity)

    # 5. Print the detailed final equation and the result
    equation_parts = [f"10^({level:.2f}/10)" for level in individual_levels]
    equation_string = " + ".join(equation_parts)
    print("The final equation to calculate the total sound level is:")
    print(f"Total Sound Level (dB) = 10 * log10({equation_string})")
    print("\nCalculating this gives:")
    print(f"Total Sound Level = {total_level_db:.2f} dB")
    
    return total_level_db

if __name__ == '__main__':
    final_answer = calculate_sound_level()
    print(f"\n<<<{final_answer:.2f}>>>")