import math

def calculate_total_sound_level():
    """
    Calculates the total sound level at a new position from multiple sources.
    """
    # Step 1 & 2: Define initial positions and the new listening position
    # My starting position is (0, 0)
    my_new_pos = {'x': 0, 'y': 25}

    sources = {
        'Dog':         {'x': -25, 'y': 0,  'L1m': 55},
        'Train':       {'x': 50,  'y': 0,  'L1m': 110},
        'Construction':{'x': 0,   'y': 75, 'L1m': 90},
        'People':      {'x': 0,   'y': -10,'L1m': 75}
    }

    print("Plan:")
    print("1. Calculate new distances from the new position (0, 25).")
    print("2. Calculate the sound level from each source at these new distances.")
    print("3. Combine the sound levels to find the total.\n")

    # This will hold the calculated levels at the new distance
    calculated_levels = {}
    total_intensity_sum = 0

    print("Calculations:")
    for name, props in sources.items():
        # Step 3: Calculate new distances
        distance = math.sqrt((props['x'] - my_new_pos['x'])**2 + (props['y'] - my_new_pos['y'])**2)

        # Step 4: Calculate sound level at the new distance
        # Formula: L2 = L1 - 20 * log10(d2 / d1). Here d1 is 1m, so log10(d2/1) = log10(d2).
        L_new = props['L1m'] - 20 * math.log10(distance)
        calculated_levels[name] = L_new

        # Prepare for Step 5: Convert to intensity and add to the sum
        total_intensity_sum += 10**(L_new / 10)

        print(f"- The {name} is {distance:.2f} meters away.")
        print(f"  The new sound level is {props['L1m']} dB - 20 * log10({distance:.2f}) = {L_new:.2f} dB\n")

    # Step 5: Combine the sound levels
    L_total = 10 * math.log10(total_intensity_sum)

    # Print the final equation with the calculated numbers
    print("Final Combined Sound Level Equation:")
    equation = f"L_total = 10 * log10(10^({calculated_levels['Dog']:.2f}/10) + 10^({calculated_levels['Train']:.2f}/10) + 10^({calculated_levels['Construction']:.2f}/10) + 10^({calculated_levels['People']:.2f}/10))"
    print(equation)
    print(f"\nTotal Sound Level = {L_total:.2f} dB")
    
    return L_total

# Run the calculation and store the final answer
final_answer = calculate_total_sound_level()
# The final result is printed above, this is just for the format requirement.
# print(f"<<<{final_answer:.2f}>>>")