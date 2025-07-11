import math

def calculate_total_sound_level():
    """
    Calculates the total sound level at a new location considering multiple sources.
    """
    # Step 1 & 2: Define initial source data and the listener's new position
    sources = {
        'dog': {'level_1m': 55, 'pos': (-25, 0)},
        'train': {'level_1m': 110, 'pos': (50, 0)},
        'construction': {'level_1m': 90, 'pos': (0, 75)},
        'people': {'level_1m': 75, 'pos': (0, -10)}
    }
    listener_new_pos = (0, 25)

    sound_intensities = []
    adjusted_levels = {}

    print("Calculating individual sound levels at the new location (25 meters towards the construction):\n")

    for name, data in sources.items():
        source_pos = data['pos']
        level_1m = data['level_1m']
        
        # Step 3: Calculate new distance to the source
        distance = math.sqrt((listener_new_pos[0] - source_pos[0])**2 + (listener_new_pos[1] - source_pos[1])**2)
        
        # Step 4: Calculate the adjusted sound level at the new distance
        # Formula: L_new = L_initial_at_1m - 20 * log10(distance)
        if distance > 0:
            adjusted_level = level_1m - 20 * math.log10(distance)
        else:
            # Handle case where distance is zero, which is unlikely but good practice
            adjusted_level = float('inf')

        adjusted_levels[name] = adjusted_level
        print(f"- {name.capitalize()}: Original Level = {level_1m} dB at 1m. New Distance = {distance:.2f}m. New Level = {adjusted_level:.2f} dB")
        
        # Convert adjusted level to relative intensity for summing
        # Formula: 10^(dB/10)
        sound_intensities.append(10**(adjusted_level / 10))

    # Step 5: Sum the intensities and convert back to a total sound level
    total_intensity = sum(sound_intensities)
    total_db = 10 * math.log10(total_intensity)

    print("\nTo find the total sound level, we combine the individual levels:")
    print(f"Total dB = 10 * log10(10^({adjusted_levels['dog']:.2f}/10) + 10^({adjusted_levels['train']:.2f}/10) + 10^({adjusted_levels['construction']:.2f}/10) + 10^({adjusted_levels['people']:.2f}/10))")
    print(f"\nFinal Combined Sound Level: {total_db:.2f} dB")
    
    return total_db

# Run the calculation and store the final answer
final_answer = calculate_total_sound_level()
# The prompt requires the answer in a specific format at the end.
# print(f"\n<<<{final_answer:.2f}>>>")
print(f'<<<{final_answer:.2f}>>>')
