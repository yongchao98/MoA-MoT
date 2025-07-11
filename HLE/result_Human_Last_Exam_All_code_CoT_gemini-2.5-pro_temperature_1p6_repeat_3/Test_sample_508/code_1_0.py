import math

def solve_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1: Define initial positions and sound levels at 1 meter.
    sources = [
        {'name': 'dog', 'pos': (-25, 0), 'L1': 55},
        {'name': 'train', 'pos': (50, 0), 'L1': 110},
        {'name': 'construction', 'pos': (0, 75), 'L1': 90},
        {'name': 'people', 'pos': (0, -10), 'L1': 75}
    ]
    reference_distance = 1.0  # r1 is 1 meter

    # Step 2: Define the observer's new position.
    my_new_pos = (0, 25)

    # Step 3 & 4: Calculate new distance and new sound level (L2) for each source.
    new_sound_levels = []
    intensity_sum_terms = []
    
    print("Calculating the sound level from each source at the new location:")
    
    for source in sources:
        # Calculate the new distance (r2)
        r2 = math.sqrt((source['pos'][0] - my_new_pos[0])**2 + (source['pos'][1] - my_new_pos[1])**2)
        
        # Calculate the new sound level (L2)
        # L2 = L1 + 20 * log10(r1 / r2)
        L2 = source['L1'] + 20 * math.log10(reference_distance / r2)
        
        new_sound_levels.append(L2)
        print(f"- The {source['name']} is {r2:.2f}m away, resulting in a sound level of {L2:.2f} dB.")

    # Step 5: Combine the sound levels.
    # L_total = 10 * log10( sum(10^(Li/10)) )
    
    total_intensity_term = sum(math.pow(10, L / 10) for L in new_sound_levels)
    
    total_db = 10 * math.log10(total_intensity_term)
    
    print("\nTo find the total sound level, we must sum the intensities, not the decibels.")
    print("The final equation is:")
    
    equation_str = "Total Sound Level (dB) = 10 * log10("
    equation_str += " + ".join([f"10^({L:.2f}/10)" for L in new_sound_levels])
    equation_str += ")"
    print(equation_str)

    print(f"\nThe calculated total sound level is: {total_db:.2f} dB")
    
    # Return final numerical answer for the system
    return total_db

if __name__ == "__main__":
    final_answer = solve_sound_level()
    # The final answer is wrapped for processing
    print(f"\n<<<{final_answer:.2f}>>>")
