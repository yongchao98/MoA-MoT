import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1 & 2: Define source locations, their sound levels at 1m, and your new position.
    sources = {
        'dog': {'pos': (-25, 0), 'Lp1': 55},
        'train': {'pos': (50, 0), 'Lp1': 110},
        'construction': {'pos': (0, 75), 'Lp1': 90},
        'people': {'pos': (0, -10), 'Lp1': 75}
    }
    my_new_position = (0, 25)
    
    new_sound_levels_db = {}
    
    print("Calculating sound levels at the new location (0, 25):\n")
    
    # Step 3 & 4: Calculate new distance and new sound level for each source.
    for name, data in sources.items():
        # Calculate the new distance from the source to the new position
        distance = math.dist(data['pos'], my_new_position)
        
        # Calculate the sound level (Lp2) at the new distance.
        # Lp2 = Lp1 - 20 * log10(r2 / r1). Here r1 = 1m.
        lp_at_new_pos = data['Lp1'] - 20 * math.log10(distance)
        new_sound_levels_db[name] = lp_at_new_pos
        
        print(f"The distance to the {name} is {distance:.2f} meters.")
        print(f"The sound level from the {name} is now {lp_at_new_pos:.2f} dB.\n")
        
    # Step 5: Calculate the total sound level by summing the intensities.
    # Total Lp = 10 * log10( Σ(10^(Lpi / 10)) )
    sum_of_intensities = sum(10**(lp / 10) for lp in new_sound_levels_db.values())
    total_db = 10 * math.log10(sum_of_intensities)

    # Format and print the final equation
    lp_dog = new_sound_levels_db['dog']
    lp_train = new_sound_levels_db['train']
    lp_construction = new_sound_levels_db['construction']
    lp_people = new_sound_levels_db['people']

    equation = (
        f"Total dB = 10 * log10(10^({lp_dog:.2f}/10) + 10^({lp_train:.2f}/10) + "
        f"10^({lp_construction:.2f}/10) + 10^({lp_people:.2f}/10))"
    )
    
    print("The total sound level is calculated by logarithmically adding the individual sound levels:")
    print(equation)
    print(f"\nThe total sound level you hear is: {total_db:.2f} dB")
    
    return total_db

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    print(f"\n<<<{final_answer:.2f}>>>")
