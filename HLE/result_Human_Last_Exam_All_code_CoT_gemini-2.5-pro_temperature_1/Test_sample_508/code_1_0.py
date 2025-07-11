import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1 & 2: Define source properties and the observer's new position
    sources = [
        {'name': 'Dog', 'level_at_1m': 55, 'pos': (-25, 0)},
        {'name': 'Train', 'level_at_1m': 110, 'pos': (50, 0)},
        {'name': 'Construction', 'level_at_1m': 90, 'pos': (0, 75)},
        {'name': 'People', 'level_at_1m': 75, 'pos': (0, -10)}
    ]
    my_new_pos = (0, 25)

    total_intensity = 0.0
    intensity_terms = []
    
    print("Calculating sound levels at the new location (0, 25):\n")

    # Step 3, 4, & 5 (part 1): Calculate distance, new dB, and intensity for each source
    for source in sources:
        source_pos = source['pos']
        
        # Step 3: Calculate new distance
        distance = math.sqrt((source_pos[0] - my_new_pos[0])**2 + (source_pos[1] - my_new_pos[1])**2)
        
        # Step 4: Calculate sound level at the new distance
        # L_new = L_at_1m - 20 * log10(distance)
        level_at_new_pos = source['level_at_1m'] - 20 * math.log10(distance)
        
        print(f"- {source['name']}:")
        print(f"  Distance = {distance:.2f} m")
        print(f"  Sound Level = {source['level_at_1m']} - 20*log10({distance:.2f}) = {level_at_new_pos:.2f} dB\n")
        
        # Step 5 (part 1): Convert dB to intensity and add to total
        intensity = 10**(level_at_new_pos / 10)
        total_intensity += intensity
        intensity_terms.append(f"10^({level_at_new_pos:.2f}/10)")

    # Step 5 (part 2): Convert total intensity back to dB
    total_db = 10 * math.log10(total_intensity)
    
    # Constructing the final equation string
    equation_str = " + ".join(intensity_terms)
    
    print("--------------------------------------------------")
    print("To find the total sound level, we sum the intensities and convert back to decibels:")
    print(f"Total dB = 10 * log10({equation_str})")
    print(f"Total dB = {total_db:.2f} dB")


if __name__ == '__main__':
    calculate_total_sound_level()
    # The final answer is derived from the calculation, so we re-calculate it here for the final output format.
    # This is a simplified re-calculation just for the final line.
    dog_dist = math.sqrt((-25-0)**2 + (0-25)**2)
    train_dist = math.sqrt((50-0)**2 + (0-25)**2)
    con_dist = math.sqrt((0-0)**2 + (75-25)**2)
    ppl_dist = math.sqrt((0-0)**2 + (-10-25)**2)
    dog_db = 55 - 20*math.log10(dog_dist)
    train_db = 110 - 20*math.log10(train_dist)
    con_db = 90 - 20*math.log10(con_dist)
    ppl_db = 75 - 20*math.log10(ppl_dist)
    total_i = 10**(dog_db/10) + 10**(train_db/10) + 10**(con_db/10) + 10**(ppl_db/10)
    final_db = 10*math.log10(total_i)
    print(f"\n<<<{final_db:.2f}>>>")