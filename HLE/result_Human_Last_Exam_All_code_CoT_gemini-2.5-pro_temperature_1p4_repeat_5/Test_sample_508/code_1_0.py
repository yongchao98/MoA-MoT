import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Initial source data: (name, level @ 1m, position (x, y))
    sources = {
        'dog': {'level': 55, 'pos': (-25, 0)},
        'train': {'level': 110, 'pos': (50, 0)},
        'construction': {'level': 90, 'pos': (0, 75)},
        'people': {'level': 75, 'pos': (0, -10)}
    }

    # 2. Listener's new position
    my_pos = (0, 25)

    # 3. Calculate attenuated sound level for each source
    attenuated_levels = {}
    intensities = {}

    print("Calculating the sound level from each source at the new location:\n")

    for name, data in sources.items():
        # Calculate the distance from the new position to the source
        distance = math.sqrt((data['pos'][0] - my_pos[0])**2 + (data['pos'][1] - my_pos[1])**2)

        # Calculate the new sound level at the listener's position
        # L2 = L1 - 20 * log10(r2 / r1). Here r1=1, so log10(r2/r1) = log10(r2)
        if distance > 0:
            level_at_listener = data['level'] - 20 * math.log10(distance)
        else:
            # Should not happen in this problem, but handle for safety
            level_at_listener = float('inf')

        attenuated_levels[name] = level_at_listener
        print(f"- {name.capitalize()}: Original Level={data['level']}dB, New Distance={distance:.2f}m, New Level={level_at_listener:.2f}dB")

    # 4. Sum the intensities and calculate the total sound level
    total_intensity = 0
    for level in attenuated_levels.values():
        total_intensity += 10**(level / 10)

    # Calculate the final total sound level in dB
    total_db = 10 * math.log10(total_intensity)

    # 5. Print the final equation and result
    l_dog = attenuated_levels['dog']
    l_train = attenuated_levels['train']
    l_construction = attenuated_levels['construction']
    l_people = attenuated_levels['people']

    print("\nTo find the total sound level, we convert each level to intensity, sum them, and convert back to decibels.")
    print("The final equation is:\n")
    print(f"Total dB = 10 * log10(10^({l_dog:.2f}/10) + 10^({l_train:.2f}/10) + 10^({l_construction:.2f}/10) + 10^({l_people:.2f}/10))")
    print("\n-------------------------------------------")
    print(f"The total calculated sound level is: {total_db:.2f} dB")
    print("-------------------------------------------\n")

    return total_db

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    print(f"<<<{final_answer:.2f}>>>")