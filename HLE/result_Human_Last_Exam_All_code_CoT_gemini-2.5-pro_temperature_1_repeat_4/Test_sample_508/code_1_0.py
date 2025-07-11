import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a new location.
    """
    # 1. Define initial positions and sound levels at 1m
    sources = {
        'dog': {'lp1': 55, 'pos': (-25, 0)},
        'train': {'lp1': 110, 'pos': (50, 0)},
        'construction': {'lp1': 90, 'pos': (0, 75)},
        'people': {'lp1': 75, 'pos': (0, -10)}
    }

    # 2. Define my new position
    my_new_pos = (0, 25)
    mx, my = my_new_pos

    new_sound_levels = {}
    total_intensity = 0

    print("Calculating sound level from each source at the new location (0, 25):\n")

    # 3 & 4. Calculate new distances and new sound levels for each source
    for name, data in sources.items():
        sx, sy = data['pos']
        lp1 = data['lp1']

        # Calculate new distance (r2)
        distance = math.sqrt((sx - mx)**2 + (sy - my)**2)

        # Calculate new sound level (Lp2)
        # Lp2 = Lp1 - 20 * log10(r2 / r1), where r1 = 1m
        lp2 = lp1 - 20 * math.log10(distance)
        new_sound_levels[name] = lp2

        # 5. Convert dB to intensity and add to the total
        intensity = 10**(lp2 / 10)
        total_intensity += intensity

        print(f"The sound level from the {name} is {lp2:.2f} dB.")

    # 5. Convert total intensity back to decibels
    total_db = 10 * math.log10(total_intensity)

    print("\nTo find the total sound level, we sum the intensities of each sound and convert back to decibels.")
    
    # Building the final equation string
    lp_dog = new_sound_levels['dog']
    lp_train = new_sound_levels['train']
    lp_construction = new_sound_levels['construction']
    lp_people = new_sound_levels['people']

    equation = (
        f"Total Sound Level = 10 * log10("
        f"10^({lp_dog:.2f}/10) + "
        f"10^({lp_train:.2f}/10) + "
        f"10^({lp_construction:.2f}/10) + "
        f"10^({lp_people:.2f}/10))"
    )

    print(f"\nFinal Equation:\n{equation}")
    print(f"\nResult:\nThe total sound level I hear is {total_db:.2f} dB.")
    
    # Return final answer for the system
    return total_db

# Execute the function
final_answer = calculate_total_sound_level()
# print(f"\n<<<{final_answer:.2f}>>>") # This is for the final answer extraction