import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial positions and sound levels (at 1 meter)
    sources = {
        'Dog':          {'coords': (-25, 0), 'db_at_1m': 55},
        'Train':        {'coords': (50, 0),  'db_at_1m': 110},
        'Construction': {'coords': (0, 75),  'db_at_1m': 90},
        'People':       {'coords': (0, -10), 'db_at_1m': 75}
    }

    # Your new position after walking 25m towards the construction site
    my_new_position = (0, 25)

    # 2. Calculate the sound level from each source at the new position
    new_db_levels = {}
    intensity_terms = []
    
    print("Calculating sound level from each source at the new location (0, 25):\n")

    for name, data in sources.items():
        source_coords = data['coords']
        db_at_1m = data['db_at_1m']

        # Calculate distance from new position to the source
        distance = math.sqrt((source_coords[0] - my_new_position[0])**2 + (source_coords[1] - my_new_position[1])**2)

        # Calculate the sound pressure level (Lp) at that distance
        # Lp_new = Lp_at_1m - 20 * log10(distance)
        lp_new = db_at_1m - 20 * math.log10(distance)
        new_db_levels[name] = lp_new
        
        # Convert dB to a value proportional to intensity for summation
        intensity = 10**(lp_new / 10)
        intensity_terms.append(intensity)

        print(f"- {name}:")
        print(f"  - Distance from (0, 25) to {source_coords} = {distance:.2f} m")
        print(f"  - Sound Level = {db_at_1m} - 20 * log10({distance:.2f}) = {lp_new:.2f} dB\n")

    # 3. Sum the intensities and convert back to a total dB level
    total_intensity = sum(intensity_terms)
    total_db = 10 * math.log10(total_intensity)

    # 4. Print the final calculation and result
    print("To find the total sound level, we convert each dB value to intensity, sum them, and convert back to dB.")
    
    lp_dog = new_db_levels['Dog']
    lp_train = new_db_levels['Train']
    lp_construction = new_db_levels['Construction']
    lp_people = new_db_levels['People']

    # The equation shows each number as requested
    final_equation = (
        f"Total dB = 10 * log10(10^({lp_dog:.2f}/10) + 10^({lp_train:.2f}/10) + "
        f"10^({lp_construction:.2f}/10) + 10^({lp_people:.2f}/10))"
    )

    print("\nFinal Equation:")
    print(final_equation)
    print(f"\nResult: {total_db:.2f} dB")
    
    return total_db

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    print(f"\n<<<{final_answer:.2f}>>>")
