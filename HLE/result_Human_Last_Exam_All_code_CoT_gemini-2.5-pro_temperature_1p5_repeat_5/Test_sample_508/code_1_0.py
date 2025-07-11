import math

def calculate_total_sound_level():
    """
    Calculates the total sound level at a new position from multiple sources.
    """
    # 1. Define initial positions and sound levels at 1 meter
    sources = {
        'dog': {'pos': (-25, 0), 'level_1m': 55},
        'train': {'pos': (50, 0), 'level_1m': 110},
        'construction': {'pos': (0, 75), 'level_1m': 90},
        'people': {'pos': (0, -10), 'level_1m': 75}
    }

    # 2. My new position after walking 25m towards construction
    my_new_pos = (0, 25)

    print("Calculating the sound level at the new location.\n")

    # 3. & 4. Calculate new distances and the sound level from each source at the new position
    new_levels = {}
    total_intensity_sum = 0
    equation_terms = []

    for name, data in sources.items():
        # Calculate new distance
        dist = math.sqrt((data['pos'][0] - my_new_pos[0])**2 + (data['pos'][1] - my_new_pos[1])**2)
        
        # Calculate sound level at the new distance
        # L_new = L_1m - 20 * log10(d_new / d_1m), where d_1m = 1
        level_at_new_pos = data['level_1m'] - 20 * math.log10(dist)
        new_levels[name] = level_at_new_pos
        
        print(f"The new distance to the {name} is {dist:.2f} meters.")
        print(f"The sound level from the {name} is {level_at_new_pos:.2f} dB.\n")
        
        # Prepare for final equation
        term_value = 10**(level_at_new_pos / 10)
        total_intensity_sum += term_value
        equation_terms.append(f"10^({level_at_new_pos:.2f}/10)")

    # 5. Calculate the total combined sound level
    total_db = 10 * math.log10(total_intensity_sum)
    
    # 6. Print the final calculation and result
    print("To find the total sound level, we use the formula:")
    print("L_total = 10 * log10(Σ(10^(Li/10)))\n")
    print("Plugging in the values for each source:")
    
    final_equation = f"L_total = 10 * log10({ ' + '.join(equation_terms) })"
    print(final_equation)
    
    print(f"\nThe total sound level you hear is: {total_db:.2f} dB")
    
    return total_db

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    print(f"\n<<<{final_answer:.2f}>>>")