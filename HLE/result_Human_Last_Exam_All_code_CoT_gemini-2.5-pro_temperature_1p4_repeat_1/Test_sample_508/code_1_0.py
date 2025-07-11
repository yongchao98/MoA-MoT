import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define source properties and your new location
    sources = [
        {'name': 'dog', 'L1': 55, 'pos': (-25, 0)},
        {'name': 'train', 'L1': 110, 'pos': (50, 0)},
        {'name': 'construction', 'L1': 90, 'pos': (0, 75)},
        {'name': 'people', 'L1': 75, 'pos': (0, -10)}
    ]
    my_pos = (0, 25)

    total_intensity = 0
    sound_levels_at_pos = {}

    print("Calculating the sound level from each source at your new location:\n")

    # 2. Calculate distance and attenuated sound level for each source
    for source in sources:
        # Calculate distance from my new position to the source
        distance = math.sqrt((source['pos'][0] - my_pos[0])**2 + (source['pos'][1] - my_pos[1])**2)

        # Calculate the sound level (in dB) from this source at my position
        # Formula: L2 = L1 - 20 * log10(d2 / d1). Here d1 is 1m, so log10(d1) is 0.
        L_at_pos = source['L1'] - 20 * math.log10(distance)
        
        sound_levels_at_pos[source['name']] = L_at_pos
        
        # 3. Convert dB to intensity and add to a running total
        intensity = 10**(L_at_pos / 10)
        total_intensity += intensity
        
        print(f"Distance to {source['name']}: {distance:.2f} meters")
        print(f"Sound level from {source['name']}: {L_at_pos:.2f} dB\n")

    # 4. Convert the total intensity back to a decibel level
    total_db = 10 * math.log10(total_intensity)
    
    # 5. Print the final calculation and result
    print("----------------------------------------------------------")
    print("To find the total sound level, we sum the sound intensities and convert back to decibels.")
    
    equation = "L_total = 10 * log10("
    equation += f"10^({sound_levels_at_pos['dog']:.2f}/10) + "  # Dog
    equation += f"10^({sound_levels_at_pos['train']:.2f}/10) + " # Train
    equation += f"10^({sound_levels_at_pos['construction']:.2f}/10) + " # Construction
    equation += f"10^({sound_levels_at_pos['people']:.2f}/10)"  # People
    equation += ")"
    
    print(f"\nFinal Equation:\n{equation}\n")
    print(f"The total combined sound level you hear is: {total_db:.2f} dB")

calculate_total_sound_level()
<<<75.12>>>