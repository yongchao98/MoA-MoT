import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a new location.
    """
    # 1. Initial setup (positions and sound levels at 1m)
    my_initial_pos = {'x': 0, 'y': 0}

    sources = {
        "Dog":          {'pos': {'x': -25, 'y': 0},  'db': 55},
        "Train":        {'pos': {'x': 50,  'y': 0},  'db': 110},
        "Construction": {'pos': {'x': 0,   'y': 75},  'db': 90},
        "People":       {'pos': {'x': 0,   'y': -10}, 'db': 75}
    }

    # 2. Determine my new position
    # Walking 25m towards the construction site (0, 75)
    my_new_pos = {'x': 0, 'y': 25}

    print("Your new position is ({}, {}).\n".format(my_new_pos['x'], my_new_pos['y']))

    # 3. & 4. Calculate new distances and resulting sound levels
    new_sound_levels = {}
    sum_of_intensities = 0
    
    print("Calculating sound level from each source at the new location:")

    for name, source in sources.items():
        # Calculate new distance
        distance = math.sqrt((source['pos']['x'] - my_new_pos['x'])**2 + (source['pos']['y'] - my_new_pos['y'])**2)
        
        # Calculate new sound level in dB
        # L_new = L_ref - 20 * log10(distance / 1)
        new_db = source['db'] - 20 * math.log10(distance)
        new_sound_levels[name] = new_db
        
        print(f"- {name}:")
        print(f"  New distance = {distance:.2f} meters")
        print(f"  Sound level = {source['db']} - 20 * log10({distance:.2f}) = {new_db:.2f} dB")
        
        # 5. Convert dB to relative intensity and add to the sum
        # I_relative = 10^(L/10)
        sum_of_intensities += 10**(new_db / 10)

    # Convert the total intensity back to dB
    # L_total = 10 * log10(Sum of I_relative)
    total_db = 10 * math.log10(sum_of_intensities)

    # Print the final calculation steps
    print("\nTo find the total sound level, we sum the relative intensities and convert back to dB.")
    print("L_total = 10 * log10(10^(L_dog/10) + 10^(L_train/10) + 10^(L_construction/10) + 10^(L_people/10))")
    print("L_total = 10 * log10(10^({:.2f}/10) + 10^({:.2f}/10) + 10^({:.2f}/10) + 10^({:.2f}/10))".format(
        new_sound_levels["Dog"], 
        new_sound_levels["Train"], 
        new_sound_levels["Construction"], 
        new_sound_levels["People"]
    ))
    print("L_total = 10 * log10({:.2f} + {:.2f} + {:.2f} + {:.2f})".format(
        10**(new_sound_levels["Dog"]/10),
        10**(new_sound_levels["Train"]/10),
        10**(new_sound_levels["Construction"]/10),
        10**(new_sound_levels["People"]/10)
    ))
    print(f"L_total = 10 * log10({sum_of_intensities:.2f})")
    print(f"L_total = {total_db:.2f} dB\n")

    print(f"The total sound level you hear is {total_db:.2f} dB.")
    return total_db

# Execute the function and capture the final answer
final_answer = calculate_total_sound_level()
print(f'<<<{final_answer:.2f}>>>')
