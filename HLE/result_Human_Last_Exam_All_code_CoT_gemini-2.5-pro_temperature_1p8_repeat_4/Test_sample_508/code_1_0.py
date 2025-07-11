import math

def calculate_total_sound_level():
    """
    Calculates the total sound level at a new location from multiple sources.
    """
    # Step 1 & 2: Define source locations and the observer's new position
    sources = [
        {'name': 'Dog', 'pos': (-25, 0), 'Lp1': 55},
        {'name': 'Train', 'pos': (50, 0), 'Lp1': 110},
        {'name': 'Construction', 'pos': (0, 75), 'Lp1': 90},
        {'name': 'People', 'pos': (0, -10), 'Lp1': 75}
    ]
    my_new_pos = (0, 25)

    new_levels_db = []
    sum_terms = []

    print("Calculating sound levels at the new location (0, 25):\n")

    for source in sources:
        # Step 3: Calculate new distance to the source
        dist_x = source['pos'][0] - my_new_pos[0]
        dist_y = source['pos'][1] - my_new_pos[1]
        new_distance = math.sqrt(dist_x**2 + dist_y**2)

        # Step 4: Calculate the new sound level (Lp2) for the source
        # Lp2 = Lp1 - 20 * log10(r2/r1), where r1=1m
        if new_distance > 0:
            new_lp = source['Lp1'] - 20 * math.log10(new_distance)
        else:
            # Handle the unlikely case of being at the source
            new_lp = float('inf')

        new_levels_db.append(new_lp)
        sum_terms.append(10**(new_lp / 10))

        print(f"Source: {source['name']}")
        print(f"  - New Distance: {new_distance:.2f} meters")
        print(f"  - Sound Level: {new_lp:.2f} dB\n")

    # Step 5: Combine the sound levels
    total_lp = 10 * math.log10(sum(sum_terms))

    # Construct and print the final equation
    equation_parts = [f"10^({level:.2f}/10)" for level in new_levels_db]
    final_equation_str = "Total Sound Level = 10 * log10(" + " + ".join(equation_parts) + ")"

    print("The final calculation is based on the formula:")
    print(final_equation_str)
    print("\n-------------------------------------------------")
    print(f"The total sound level you hear is: {total_lp:.2f} dB")
    print("-------------------------------------------------")
    
    return total_lp

if __name__ == '__main__':
    final_answer = calculate_total_sound_level()
    print(f'<<<{final_answer:.2f}>>>')
