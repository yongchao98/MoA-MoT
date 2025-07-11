import math

def calculate_total_sound_level():
    """
    Calculates the total sound level at a new location based on multiple sound sources.
    """
    # 1. Define initial positions and sound source data
    # My initial position is (0,0). I walk 25m towards the construction site.
    my_new_position = (0, 25)

    sources = [
        {'name': 'Dog', 'level_at_1m': 55, 'position': (-25, 0)},
        {'name': 'Train', 'level_at_1m': 110, 'position': (50, 0)},
        {'name': 'Construction', 'level_at_1m': 90, 'position': (0, 75)},
        {'name': 'People', 'level_at_1m': 75, 'position': (0, -10)}
    ]

    new_levels_db = []
    sum_of_intensities = 0

    print("Calculating the sound level from each source at the new location (0, 25):\n")

    # 2. & 3. Calculate new distance and attenuated sound level for each source
    for source in sources:
        x1, y1 = source['position']
        x2, y2 = my_new_position
        
        # Calculate distance from the new position to the source
        distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        # Calculate the sound level (L2) at the new distance (d2)
        # L2 = L1 - 20 * log10(d2 / d1), where d1 = 1 meter
        level_at_distance = source['level_at_1m'] - 20 * math.log10(distance)
        new_levels_db.append(level_at_distance)
        
        print(f"Source: {source['name']}")
        print(f"  - Distance from new location: {distance:.2f} meters")
        print(f"  - Sound level at new location: {source['level_at_1m']} dB - 20*log10({distance:.2f}) = {level_at_distance:.2f} dB\n")
        
        # 4. Convert dB to intensity and add to the sum
        sum_of_intensities += 10**(level_at_distance / 10)
        
    # 5. Convert total intensity back to a total dB level
    total_db = 10 * math.log10(sum_of_intensities)
    
    print("--------------------------------------------------")
    print("To find the total sound level, we sum the intensities (not the decibels).")
    print("The formula is: L_total = 10 * log10( Σ(10^(Li/10)) )")
    print("\nFinal Calculation:")
    equation = f"Total dB = 10 * log10(10^({new_levels_db[0]:.2f}/10) + 10^({new_levels_db[1]:.2f}/10) + 10^({new_levels_db[2]:.2f}/10) + 10^({new_levels_db[3]:.2f}/10))"
    print(equation)
    
    print(f"\nTotal Sound Level: {total_db:.2f} dB")
    
    return total_db

# Run the calculation and store the result
final_answer = calculate_total_sound_level()

# Print the final answer in the specified format
print(f"<<<{final_answer:.2f}>>>")