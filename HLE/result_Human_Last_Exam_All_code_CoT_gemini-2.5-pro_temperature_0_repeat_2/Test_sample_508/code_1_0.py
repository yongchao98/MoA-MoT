import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1 & 2: Define source properties and the observer's new position
    # Sources are defined as a dictionary:
    # key: name
    # value: dictionary with sound level at 1m (lp1) and position (x, y)
    sources = {
        "dog": {"lp1": 55, "pos": (-25, 0)},
        "train": {"lp1": 110, "pos": (50, 0)},
        "construction": {"lp1": 90, "pos": (0, 75)},
        "people": {"lp1": 75, "pos": (0, -10)}
    }
    my_new_pos = (0, 25)

    new_sound_levels = {}
    sound_intensities = []

    # Step 3 & 4: Calculate new distance and resulting sound level for each source
    print("Calculating the sound level from each source at the new location:")
    for name, data in sources.items():
        source_pos = data["pos"]
        lp1 = data["lp1"]
        
        # Calculate the new distance from the observer to the source
        distance = math.sqrt((source_pos[0] - my_new_pos[0])**2 + (source_pos[1] - my_new_pos[1])**2)
        
        # Calculate the sound level at the new distance
        # Formula: Lp2 = Lp1 - 20 * log10(r2/r1), where r1 is 1m
        lp2 = lp1 - 20 * math.log10(distance)
        new_sound_levels[name] = lp2
        
        print(f"- The sound level from the {name} (originally {lp1} dB at 1m) at a new distance of {distance:.2f}m is {lp2:.2f} dB.")

    # Step 5: Sum the sound levels
    # Convert each dB level to intensity, sum them, and convert back to dB
    # Formula: L_total = 10 * log10(sum(10^(Li/10)))
    for level in new_sound_levels.values():
        sound_intensities.append(10**(level / 10))
        
    total_intensity = sum(sound_intensities)
    total_sound_level = 10 * math.log10(total_intensity)

    # Unpack the calculated levels for the final equation string
    lp_dog = new_sound_levels['dog']
    lp_train = new_sound_levels['train']
    lp_construction = new_sound_levels['construction']
    lp_people = new_sound_levels['people']

    print("\nThe total sound level is calculated by summing the sound intensities:")
    print(f"Total Sound Level = 10 * log10(10^({lp_dog:.2f}/10) + 10^({lp_train:.2f}/10) + 10^({lp_construction:.2f}/10) + 10^({lp_people:.2f}/10))")
    print(f"\nThe final calculated total sound level is: {total_sound_level:.2f} dB")
    
    return total_sound_level

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    # The final answer is wrapped for programmatic retrieval.
    print(f"\n<<<{final_answer:.2f}>>>")
