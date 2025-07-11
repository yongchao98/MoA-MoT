import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1 & 2: Define source positions and the observer's new position
    # Format: {'name': 'Source Name', 'L1': dB at 1m, 'pos': (x, y)}
    sources = [
        {'name': 'dog', 'L1': 55, 'pos': (-25, 0)},
        {'name': 'train', 'L1': 110, 'pos': (50, 0)},
        {'name': 'construction', 'L1': 90, 'pos': (0, 75)},
        {'name': 'people', 'L1': 75, 'pos': (0, -10)}
    ]
    my_pos = (0, 25)

    total_intensity = 0
    calculated_levels = {}

    print("Calculating the sound level from each source at your new location:")
    # Step 3, 4 & 5 (part 1): Calculate distance, new level, and intensity for each source
    for source in sources:
        # Step 3: Calculate distance
        dist = math.sqrt((source['pos'][0] - my_pos[0])**2 + (source['pos'][1] - my_pos[1])**2)

        # Step 4: Calculate sound level (L2) at the new distance
        # L2 = L1 - 20 * log10(d2/d1), where d1 is 1 meter.
        L2 = source['L1'] - 20 * math.log10(dist)
        calculated_levels[source['name']] = L2
        print(f"- {source['name'].capitalize()}: {L2:.2f} dB (at a distance of {dist:.2f} meters)")

        # Step 5 (part 1): Convert dB to intensity and add to total
        intensity = 10**(L2 / 10)
        total_intensity += intensity

    # Step 6: Convert total intensity back to a total sound level in dB
    total_db = 10 * math.log10(total_intensity)

    # Output the final equation with the calculated numbers
    dog_db = calculated_levels['dog']
    train_db = calculated_levels['train']
    construction_db = calculated_levels['construction']
    people_db = calculated_levels['people']

    print("\nTo find the total sound level, we sum the intensities (not the dB values directly).")
    print("The final equation is:")
    print(f"Total dB = 10 * log10(10^({dog_db:.2f}/10) + 10^({train_db:.2f}/10) + 10^({construction_db:.2f}/10) + 10^({people_db:.2f}/10))")
    print(f"\nThe total sound level you hear is {total_db:.2f} dB.")

if __name__ == '__main__':
    calculate_total_sound_level()
<<<75.12>>>