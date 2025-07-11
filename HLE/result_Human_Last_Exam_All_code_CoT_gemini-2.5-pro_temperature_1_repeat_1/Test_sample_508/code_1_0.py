import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1 & 2: Define initial source properties and the observer's new position.
    # Format: { 'name': [L1 (dB at 1m), x_pos, y_pos] }
    sources = {
        'Dog': [55, -25, 0],
        'Train': [110, 50, 0],
        'Construction': [90, 0, 75],
        'People': [75, 0, -10]
    }
    observer_pos = (0, 25)

    print("Observer's new position: (0, 25)\n")

    # Step 3 & 4: Calculate new distance and new sound level for each source.
    new_levels_db = {}
    print("Calculating the new sound level from each source at the observer's location:")
    for name, properties in sources.items():
        l1 = properties[0]
        source_pos = (properties[1], properties[2])

        # Calculate new distance (r2)
        r2 = math.sqrt((source_pos[0] - observer_pos[0])**2 + (source_pos[1] - observer_pos[1])**2)

        # Calculate new sound level (L2)
        # L2 = L1 - 20 * log10(r2 / r1), where r1 = 1m
        l2 = l1 - 20 * math.log10(r2)
        new_levels_db[name] = l2

        print(f"- {name}:")
        print(f"  - New distance: {r2:.2f} meters")
        print(f"  - New sound level = {l1} dB - 20 * log10({r2:.2f}) = {l2:.2f} dB")

    # Step 5: Combine the sound levels.
    # L_total = 10 * log10( sum(10^(Li/10)) )
    sum_of_intensities = 0
    for l2 in new_levels_db.values():
        sum_of_intensities += 10**(l2 / 10)

    total_level_db = 10 * math.log10(sum_of_intensities)

    # Print the final calculation equation
    print("\nCombining the sound levels:")
    equation_parts = []
    for name, l2 in new_levels_db.items():
        equation_parts.append(f"10^({l2:.2f}/10)")

    final_equation = "Total Sound Level = 10 * log10(" + " + ".join(equation_parts) + ")"
    print(final_equation)

    print(f"\nFinal calculated total sound level: {total_level_db:.2f} dB")
    
    return total_level_db

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    print(f"\n<<<{final_answer:.2f}>>>")
