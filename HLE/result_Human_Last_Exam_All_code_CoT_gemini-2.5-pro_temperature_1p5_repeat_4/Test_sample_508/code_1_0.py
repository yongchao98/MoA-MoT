import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a new location.
    """
    # 1. Define initial positions and sound levels at 1 meter.
    # My initial position is the origin (0, 0).
    # Front is positive y, Right is positive x.
    sources = {
        'dog': {'pos': (-25, 0), 'lp1': 55},
        'train': {'pos': (50, 0), 'lp1': 110},
        'construction': {'pos': (0, 75), 'lp1': 90},
        'people': {'pos': (0, -10), 'lp1': 75}
    }

    # 2. Define the new position after walking 25m towards the construction site.
    my_new_pos = (0, 25)

    new_sound_levels = {}
    intensities = []

    print("Calculating sound levels at the new position (0, 25):\n")

    # 3. & 4. Calculate new distances and the corresponding sound level for each source.
    for name, data in sources.items():
        # Calculate distance from my new position to the source
        dist = math.sqrt((data['pos'][0] - my_new_pos[0])**2 + (data['pos'][1] - my_new_pos[1])**2)
        
        # Calculate the sound level (Lp) at the new distance.
        # Lp2 = Lp1 - 20 * log10(r2 / r1), where r1=1m.
        # This simplifies to Lp2 = Lp1 - 20 * log10(r2).
        if dist > 0:
            lp_new = data['lp1'] - 20 * math.log10(dist)
        else:
            # Cannot be at distance 0, but as a fallback, avoid log(0) error.
            # In a real-world scenario, you wouldn't be at the exact point source.
            lp_new = float('inf') 
        
        new_sound_levels[name] = lp_new
        print(f"Distance to {name}: {dist:.2f} meters, New Sound Level: {lp_new:.2f} dB")

    # 5. Combine the sound levels.
    # L_total = 10 * log10(sum of 10^(Li/10))
    sum_of_intensities = 0
    for level in new_sound_levels.values():
        sum_of_intensities += 10**(level / 10)
    
    total_db = 10 * math.log10(sum_of_intensities)

    # Print the final equation with calculated values
    l_dog = new_sound_levels['dog']
    l_train = new_sound_levels['train']
    l_construction = new_sound_levels['construction']
    l_people = new_sound_levels['people']

    print("\n------------------------------------------------------------")
    print("The total sound level is the sum of the sound intensities from each source.")
    print("The formula is: L_total = 10 * log10(10^(L1/10) + 10^(L2/10) + 10^(L3/10) + 10^(L4/10))")
    print("\nFinal calculation:")
    final_equation = f"Total Level = 10 * log10(10^({l_dog:.2f}/10) + 10^({l_train:.2f}/10) + 10^({l_construction:.2f}/10) + 10^({l_people:.2f}/10))"
    print(final_equation)
    print(f"\nThe total sound level you hear is: {total_db:.2f} dB")
    print("------------------------------------------------------------")

    return total_db

# Run the calculation and store the result
final_sound_level = calculate_total_sound_level()

# Final answer in the required format
print(f"\n<<<{final_sound_level:.2f}>>>")
