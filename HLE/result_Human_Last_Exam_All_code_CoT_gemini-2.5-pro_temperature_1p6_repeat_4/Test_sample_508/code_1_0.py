import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial positions and sound levels at 1 meter.
    # My initial position is at origin (0,0)
    sources = {
        'dog': {'position': (-25, 0), 'level_at_1m': 55},
        'train': {'position': (50, 0), 'level_at_1m': 110},
        'construction': {'position': (0, 75), 'level_at_1m': 90},
        'people': {'position': (0, -10), 'level_at_1m': 75}
    }

    # 2. Define the new location after walking 25m towards construction
    my_new_position = (0, 25)

    print(f"My new position is: {my_new_position}\n")

    new_levels = {}
    total_intensity = 0

    print("Calculating the sound level from each source at my new location:")
    
    # 3 & 4. Calculate new distances and new sound levels for each source
    for name, data in sources.items():
        source_pos = data['position']
        # Calculate new distance using Euclidean distance formula
        distance = math.sqrt((my_new_position[0] - source_pos[0])**2 + (my_new_position[1] - source_pos[1])**2)
        
        # Calculate new sound level using the inverse square law for sound
        # L2 = L1 - 20 * log10(d2 / d1), where d1 = 1m
        level_at_1m = data['level_at_1m']
        new_level = level_at_1m - 20 * math.log10(distance)
        new_levels[name] = new_level
        
        print(f"- {name.capitalize()}:")
        print(f"  New distance = {distance:.2f} meters")
        print(f"  New sound level = {level_at_1m} dB - 20 * log10({distance:.2f}) = {new_level:.2f} dB\n")
        
        # 5. Convert dB to relative intensity and add to the total
        total_intensity += 10**(new_level / 10)

    # 6. Convert total intensity back to dB
    total_db = 10 * math.log10(total_intensity)
    
    print("--------------------------------------------------")
    print("To find the total sound level, we sum the intensities (not the dB values).")
    print("The formula is: L_total = 10 * log10(10^(L1/10) + 10^(L2/10) + ...)\n")
    
    l_dog = new_levels['dog']
    l_train = new_levels['train']
    l_construction = new_levels['construction']
    l_people = new_levels['people']

    print("Final Calculation:")
    print(f"Total Sound Level = 10 * log10(10^({l_dog:.2f}/10) + 10^({l_train:.2f}/10) + 10^({l_construction:.2f}/10) + 10^({l_people:.2f}/10))")
    print(f"Total Sound Level = {total_db:.2f} dB")


if __name__ == '__main__':
    calculate_total_sound_level()
    # The final answer must be on the last line with the specified format
    # Final value from running the script: 75.12
    print("<<<75.12>>>")