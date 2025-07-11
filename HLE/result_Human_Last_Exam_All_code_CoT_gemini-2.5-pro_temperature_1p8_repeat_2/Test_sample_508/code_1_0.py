import math

def solve_sound_level():
    """
    Calculates the total sound level from multiple sources at a new location.
    """
    # 1. Define sources: name, reference sound level (dB at 1m), and position (x, y)
    sources = {
        'Dog': {'L1': 55, 'pos': (-25, 0)},
        'Train': {'L1': 110, 'pos': (50, 0)},
        'Construction': {'L1': 90, 'pos': (0, 75)},
        'People': {'L1': 75, 'pos': (0, -10)}
    }

    # 2. Define your new position
    my_pos = (0, 25)

    print("Calculating the sound level at the new position (0, 25).\n")

    # Dictionary to store the calculated sound level for each source at the listener's position
    attenuated_levels = {}
    
    # List to store the terms for the final summation equation
    summation_terms = []

    for name, data in sources.items():
        source_pos = data['pos']
        l1 = data['L1']

        # 3. Calculate the distance from the new position to the source
        distance = math.sqrt((source_pos[0] - my_pos[0])**2 + (source_pos[1] - my_pos[1])**2)

        # 4. Calculate the attenuated sound level (L2) at the new distance
        # L2 = L1 - 20 * log10(d2 / d1), where d1 = 1m
        if distance < 1:
            # Avoid log(d) for d<1 which would increase the sound level.
            # Practically, at distances closer than the reference, the level would be higher,
            # but for this problem, we assume the listener is not that close.
            # If d=0, it's a singularity. We'll cap at d=1 for this calculation.
             distance = 1

        l2 = l1 - 20 * math.log10(distance)
        attenuated_levels[name] = l2
        
        summation_terms.append(f"10^({l2:.2f}/10)")

        print(f"Source: {name}")
        print(f"  - Distance from new location: {distance:.2f} meters")
        print(f"  - Sound level at new location: {l1} - 20 * log10({distance:.2f}) = {l2:.2f} dB\n")

    # 5. Combine the sound levels
    # Convert each dB level to relative intensity (10^(L/10)), sum them up, and convert back to dB.
    total_intensity = sum(10**(level / 10) for level in attenuated_levels.values())
    total_db = 10 * math.log10(total_intensity)

    # Print the final calculation step
    print("Combining the sound levels:")
    equation = " + ".join(summation_terms)
    print(f"Total Level (dB) = 10 * log10({equation})")
    
    final_equation_values = " + ".join([f"{10**(l/10):.2f}" for l in attenuated_levels.values()])
    print(f"Total Level (dB) = 10 * log10({final_equation_values})")
    print(f"Total Level (dB) = 10 * log10({total_intensity:.2f})")
    
    print("\n--- Final Answer ---")
    print(f"The total sound level you hear is: {total_db:.2f} dB")
    
    # For the final answer tag format
    return total_db

# Run the solver
final_answer = solve_sound_level()
# The final numerical result needs to be printed for the answer format below
# print(f'<<<{final_answer:.2f}>>>')