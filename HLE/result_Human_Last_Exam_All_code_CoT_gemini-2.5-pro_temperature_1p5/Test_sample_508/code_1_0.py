import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial parameters and positions
    my_position = (0, 25)  # Moved 25m towards the construction site

    sources = {
        'Dog':      {'level_at_1m': 55, 'pos': (-25, 0)},
        'Train':    {'level_at_1m': 110, 'pos': (50, 0)},
        'Construction': {'level_at_1m': 90, 'pos': (0, 75)},
        'People':   {'level_at_1m': 75, 'pos': (0, -10)}
    }

    print("Calculating the sound level from each source at your new location (0, 25):\n")

    # 2. Calculate the new sound level for each source
    new_levels_db = {}
    intensities = {}
    total_intensity = 0.0

    for name, data in sources.items():
        # Calculate distance from new position to the source
        dist = math.sqrt((data['pos'][0] - my_position[0])**2 + (data['pos'][1] - my_position[1])**2)
        
        # Calculate new sound level (L2 = L1 - 20*log10(d2)) as d1 is 1m
        # Clamp to 0 if distance is less than or equal to 1m to avoid log(d)<=0 issues
        if dist > 1:
            new_level = data['level_at_1m'] - 20 * math.log10(dist)
        else: # If distance is <= 1m, the sound doesn't attenuate by this formula
             new_level = data['level_at_1m']
        
        new_levels_db[name] = new_level
        
        # Convert dB to relative intensity (I = 10^(L/10))
        intensity = 10**(new_level / 10)
        intensities[name] = intensity
        
        # Add to total intensity
        total_intensity += intensity
        
        print(f"- {name}:")
        print(f"  Distance = {dist:.2f} meters")
        print(f"  Sound Level = {data['level_at_1m']} dB - 20*log10({dist:.2f}) = {new_level:.2f} dB\n")

    # 3. Convert total intensity back to a total sound level
    total_sound_level = 10 * math.log10(total_intensity)

    # 4. Print the final calculation and result
    print("To find the total sound level, we convert each dB level to intensity, sum them up, and convert back to dB.")
    print("Total Sound Level = 10 * log10(10^(L_dog/10) + 10^(L_train/10) + 10^(L_construction/10) + 10^(L_people/10))")
    
    final_eq = "Total Sound Level = 10 * log10("
    final_eq += f"10^({new_levels_db['Dog']:.2f}/10)"
    final_eq += f" + 10^({new_levels_db['Train']:.2f}/10)"
    final_eq += f" + 10^({new_levels_db['Construction']:.2f}/10)"
    final_eq += f" + 10^({new_levels_db['People']:.2f}/10)"
    final_eq += f") = {total_sound_level:.2f} dB"
    print(final_eq)

    print(f"\nThe total sound level you hear is: {total_sound_level:.2f} dB")
    
    return total_sound_level

# Execute the calculation
final_answer = calculate_total_sound_level()
# The final result is wrapped for programmatic access.
# print(f"\n<<<{final_answer:.2f}>>>")
# The user wants me to not output the <<< part in a separate print statement
# but directly as part of the response text. So I'll modify the final part.

if __name__ == '__main__':
    pass # No need to call again, already called.