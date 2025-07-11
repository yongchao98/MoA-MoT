import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial sources and my position
    sources = {
        'dog': {'L1': 55, 'pos': (-25, 0)},
        'train': {'L1': 110, 'pos': (50, 0)},
        'construction': {'L1': 90, 'pos': (0, 75)},
        'people': {'L1': 75, 'pos': (0, -10)}
    }
    my_new_pos = (0, 25)
    r1 = 1.0  # Reference distance is 1 meter

    print("Calculating the total sound level at the new position (0, 25).\n")

    total_intensity = 0
    calculated_levels = {}

    # 2-4. Calculate new distance, new sound level, and intensity for each source
    for name, data in sources.items():
        # Calculate new distance (r2)
        r2 = math.sqrt((my_new_pos[0] - data['pos'][0])**2 + (my_new_pos[1] - data['pos'][1])**2)
        
        # Calculate new sound level (L2)
        # L2 = L1 - 20 * log10(r2 / r1). Since r1=1, this simplifies to L1 - 20 * log10(r2)
        l2 = data['L1'] - 20 * math.log10(r2 / r1)
        calculated_levels[name] = l2
        
        # Convert L2 to intensity and add to total
        intensity = 10**(l2 / 10)
        total_intensity += intensity
        
        print(f"Source: {name.capitalize()}")
        print(f"  - Distance to new position: {r2:.2f} meters")
        print(f"  - Sound level at new position: {l2:.2f} dB")
        print("-" * 20)

    # Convert total intensity back to a total sound level
    total_sound_level = 10 * math.log10(total_intensity)

    print("\nTo find the total sound level, we convert each source's dB level to intensity,")
    print("sum the intensities, and convert the result back to dB.")
    
    # Building the final equation string with each calculated level
    dog_l = calculated_levels['dog']
    train_l = calculated_levels['train']
    construction_l = calculated_levels['construction']
    people_l = calculated_levels['people']

    equation = (
        f"Total dB = 10 * log10(10^({dog_l:.2f}/10) + 10^({train_l:.2f}/10) + "
        f"10^({construction_l:.2f}/10) + 10^({people_l:.2f}/10))"
    )

    print("\nFinal Equation:")
    print(equation)
    
    print(f"\nThe total sound level you hear is: {total_sound_level:.2f} dB")
    return total_sound_level

# Execute the function and capture the final numerical answer
final_answer = calculate_total_sound_level()
print(f"\n<<<{final_answer:.2f}>>>")
