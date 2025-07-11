import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a new location.
    """
    # Step 1 & 2: Define initial source positions, reference sound levels, and the new observer position.
    # Format: {'name': str, 'Lp_ref': dB, 'pos': (x, y)}
    # Lp_ref is the sound level measured at 1 meter.
    sources = [
        {'name': 'Dog', 'Lp_ref': 55, 'pos': (-25, 0)},
        {'name': 'Train', 'Lp_ref': 110, 'pos': (50, 0)},
        {'name': 'Construction', 'Lp_ref': 90, 'pos': (0, 75)},
        {'name': 'People', 'Lp_ref': 75, 'pos': (0, -10)}
    ]
    my_new_pos = (0, 25)

    new_sound_levels = []
    sum_of_powers = 0

    print("Calculating individual sound levels at the new position:")
    print("-" * 50)
    
    # Step 3 & 4: Calculate new distances and the new sound level for each source.
    for source in sources:
        x1, y1 = source['pos']
        x2, y2 = my_new_pos
        
        distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
        
        # Formula: Lp2 = Lp1 - 20 * log10(d2/d1). Here d1 is 1m.
        Lp_new = source['Lp_ref'] - 20 * math.log10(distance)
        
        new_sound_levels.append(Lp_new)
        print(f"Sound level from {source['name']:<12} at {distance:>5.2f}m away: {Lp_new:>5.2f} dB")

    # Step 5: Sum the sound levels.
    for lp in new_sound_levels:
        sum_of_powers += 10**(lp / 10)

    total_Lp = 10 * math.log10(sum_of_powers)

    print("-" * 50)
    print("To find the total sound level, we use the formula:")
    print("L_total = 10 * log10(10^(L1/10) + 10^(L2/10) + ... + 10^(Ln/10))\n")
    
    equation_parts = [f"10^({lp:.2f}/10)" for lp in new_sound_levels]
    print(f"L_total = 10 * log10({ ' + '.join(equation_parts) })")
    print(f"\nTotal Sound Level at the new location: {total_Lp:.2f} dB")
    
    return total_Lp

# Execute the calculation and capture the final answer.
final_answer = calculate_total_sound_level()
print(f'<<<{final_answer:.2f}>>>')