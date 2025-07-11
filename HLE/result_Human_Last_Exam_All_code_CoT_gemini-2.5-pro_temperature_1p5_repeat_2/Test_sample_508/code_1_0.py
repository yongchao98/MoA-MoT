import math

def calculate_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Initial positions and sound levels (L1 at r1=1m)
    sources = {
        'Dog':          {'L1': 55,  'pos': (-25, 0)},
        'Train':        {'L1': 110, 'pos': (50, 0)},
        'Construction': {'L1': 90,  'pos': (0, 75)},
        'People':       {'L1': 75,  'pos': (0, -10)}
    }

    # 2. Your new position after walking 25m towards the construction site
    my_pos = (0, 25)

    # Dictionary to store calculated values for each source
    results = {}
    
    # 3. & 4. Calculate new distances (r2) and new sound levels (L2)
    print("Calculating the sound level (dB) from each source at the new location:\n")
    for name, data in sources.items():
        # Calculate new distance r2
        dist_x = data['pos'][0] - my_pos[0]
        dist_y = data['pos'][1] - my_pos[1]
        r2 = math.sqrt(dist_x**2 + dist_y**2)

        # Calculate new sound level L2 = L1 - 20 * log10(r2 / r1), where r1=1
        if r2 < 1: # Avoid log(negative) if standing on top of the source
            r2 = 1 
        L2 = data['L1'] - 20 * math.log10(r2)
        
        # Store results
        results[name] = {'r2': r2, 'L2': L2}
        print(f"- {name}: New distance is {r2:.2f} m, new sound level is {L2:.2f} dB.")

    # 5. Combine the sound levels
    print("\nTo find the total sound level, we convert each dB value to intensity, sum them, and convert back to dB.")
    print("Formula: Total dB = 10 * log10( 10^(L_dog/10) + 10^(L_train/10) + ... )\n")

    # Calculate relative intensities
    total_intensity = 0
    intensity_eq_parts = []
    db_eq_parts = []
    
    for name, data in results.items():
        intensity = 10**(data['L2'] / 10)
        results[name]['intensity'] = intensity
        intensity_eq_parts.append(f"{intensity:.1f}")
        db_eq_parts.append(f"10^({data['L2']:.2f}/10)")

    # Sum the intensities
    total_intensity = sum(data['intensity'] for data in results.values())
    
    # Calculate total dB level
    total_db = 10 * math.log10(total_intensity)

    # 6. Print the final calculation steps
    print("The final calculation is:")
    print(f"Total dB = 10 * log10( {' + '.join(db_eq_parts)} )")
    print(f"Total dB = 10 * log10( {' + '.join(intensity_eq_parts)} )")
    print(f"Total dB = 10 * log10( {total_intensity:.1f} )")
    print(f"Total dB = {10 * math.log10(total_intensity):.2f}\n")
    
    print(f"The total sound level you hear is {total_db:.2f} dB.")
    
    return total_db

if __name__ == "__main__":
    final_answer = calculate_sound_level()
    print(f"\n<<<{final_answer:.2f}>>>")
