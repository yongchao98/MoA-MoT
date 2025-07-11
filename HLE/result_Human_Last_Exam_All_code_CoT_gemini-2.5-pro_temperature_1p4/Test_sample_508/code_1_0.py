import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial positions and sound source data
    # My new position after walking 25m towards the construction site
    my_new_position = (0, 25)

    sources = {
        'dog': {'L_ref': 55, 'pos': (-25, 0)},
        'train': {'L_ref': 110, 'pos': (50, 0)},
        'construction': {'L_ref': 90, 'pos': (0, 75)},
        'people': {'L_ref': 75, 'pos': (0, -10)}
    }

    new_levels_db = {}
    total_intensity_sum = 0

    print("Calculating the sound level from each source at the new location (0, 25):\n")

    # 2. & 3. Calculate new distances and sound levels for each source
    for name, data in sources.items():
        source_pos = data['pos']
        L_ref = data['L_ref']

        # Calculate distance from my new position to the source
        distance = math.sqrt((source_pos[0] - my_new_position[0])**2 + (source_pos[1] - my_new_position[1])**2)

        # Calculate the new sound level (in dB) at my position
        # L_new = L_ref - 20 * log10(distance)
        L_new = L_ref - 20 * math.log10(distance)
        new_levels_db[name] = L_new

        # 4. Convert dB to intensity term and add to sum
        total_intensity_sum += 10**(L_new / 10)

        print(f"Distance to {name}: {distance:.2f} meters")
        print(f"Sound level from {name}: {L_new:.2f} dB\n")

    # 5. Combine the sound levels
    # L_total = 10 * log10(sum of all 10^(L_i/10))
    total_db = 10 * math.log10(total_intensity_sum)
    
    # Print the final equation with the calculated numbers
    dog_db = new_levels_db['dog']
    train_db = new_levels_db['train']
    construction_db = new_levels_db['construction']
    people_db = new_levels_db['people']
    
    print("--------------------------------------------------")
    print("The total sound level is calculated by combining the intensities of each sound source.")
    print("The formula is: L_total = 10 * log10(10^(L1/10) + 10^(L2/10) + ...)\n")
    print("Plugging in the numbers for each source:")
    print(f"L_total = 10 * log10(10^({dog_db:.2f}/10) + 10^({train_db:.2f}/10) + 10^({construction_db:.2f}/10) + 10^({people_db:.2f}/10))")
    
    print("\n--------------------------------------------------")
    print(f"The total sound level you hear is: {total_db:.2f} dB")
    
    return total_db

# Execute the function
total_sound_level = calculate_total_sound_level()
# The final answer in the required format
# print(f"<<<{total_sound_level:.2f}>>>")
# The problem asks for one final answer, let's round to one decimal place as is common for dB.
print(f"\n<<<{round(total_sound_level, 1)}>>>")
