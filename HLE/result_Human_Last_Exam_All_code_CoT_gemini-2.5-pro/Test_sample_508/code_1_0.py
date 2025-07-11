import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial positions and sound source properties
    # My new position: (0, 25)
    my_pos = (0, 25)

    # Sound sources: name, level at 1m (dB), initial position
    sources = [
        {'name': 'Dog', 'L1': 55, 'pos': (-25, 0)},
        {'name': 'Train', 'L1': 110, 'pos': (50, 0)},
        {'name': 'Construction', 'L1': 90, 'pos': (0, 75)},
        {'name': 'People', 'L1': 75, 'pos': (0, -10)}
    ]

    individual_levels = []
    total_intensity = 0
    equation_parts = []

    print("Calculating the sound level for each source at your new location:\n")

    # 2. Calculate new distance and sound level for each source
    for source in sources:
        # Calculate distance from my new position to the source
        distance = math.sqrt((source['pos'][0] - my_pos[0])**2 + (source['pos'][1] - my_pos[1])**2)

        # Calculate sound level (L2) at the new distance
        # Formula: L2 = L1 - 20 * log10(distance)
        level_at_listener = source['L1'] - 20 * math.log10(distance)
        individual_levels.append(level_at_listener)

        print(f"{source['name']}:")
        print(f"  - Distance to new location: {distance:.2f} meters")
        print(f"  - Sound level at new location: {level_at_listener:.2f} dB\n")

        # 3. Convert level to intensity and add to total
        # Formula: I = 10^(L/10)
        intensity = 10**(level_at_listener / 10)
        total_intensity += intensity

        # Prepare parts for the final equation string
        equation_parts.append(f"10^({level_at_listener:.2f}/10)")

    # 4. Calculate total sound level from total intensity
    # Formula: L_total = 10 * log10(I_total)
    total_db = 10 * math.log10(total_intensity)

    # 5. Print the final results and the formula used
    print("----------------------------------------------------------")
    print("To find the combined sound level, we sum the sound intensities, not the decibels.")
    print("The formula is: Total dB = 10 * log10( Σ(10^(Li/10)) )")
    print("\nApplying this to our calculated levels:")
    
    equation_string = " + ".join(equation_parts)
    print(f"Total dB = 10 * log10({equation_string})")

    print("\n----------------------------------------------------------")
    print(f"The total sound level you hear is: {total_db:.2f} dB")
    print("----------------------------------------------------------")
    
    return total_db

# Run the calculation and store the final answer
final_answer = calculate_total_sound_level()
# The final answer in the required format
# print(f"<<<{final_answer:.2f}>>>")