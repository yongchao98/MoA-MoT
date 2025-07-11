import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a new location.
    """
    # Step 1 & 2: Define initial positions and the new position
    # My initial position is (0,0). My new position is (0, 25).
    my_pos = (0, 25)

    sources = {
        "Dog":          {"pos": (-25, 0), "db_at_1m": 55},
        "Train":        {"pos": (50, 0),  "db_at_1m": 110},
        "Construction": {"pos": (0, 75),  "db_at_1m": 90},
        "People":       {"pos": (0, -10), "db_at_1m": 75}
    }

    new_db_levels = {}
    intensities = {}
    total_intensity = 0

    print("Calculating sound levels at the new position (0, 25):\n")

    # Step 3 & 4: Calculate new distance and new dB level for each source
    for name, data in sources.items():
        # Calculate new distance
        dist = math.sqrt((data["pos"][0] - my_pos[0])**2 + (data["pos"][1] - my_pos[1])**2)
        
        # Calculate new sound level (L_new = L_at_1m - 20 * log10(distance))
        # We handle distance < 1m case, though not present in this problem.
        if dist < 1:
            # This formula is for d >= 1. If d < 1, sound would increase.
            # L_new = L_at_1m + 20 * log10(1/dist)
            # For simplicity, we assume distance is always >= 1m from the source center.
            # In this problem, all distances are > 1.
            pass
        
        db_new = data["db_at_1m"] - 20 * math.log10(dist)
        new_db_levels[name] = db_new
        
        print(f"- {name}:")
        print(f"  New distance = {dist:.2f} meters")
        print(f"  Sound level = {data['db_at_1m']} - 20*log10({dist:.2f}) = {db_new:.2f} dB")

    # Step 5: Convert to intensities, sum them, and convert back to total dB
    print("\nCombining the sound levels:")
    
    # Build the equation string
    equation_parts = []
    for name, db in new_db_levels.items():
        intensity = 10**(db / 10)
        intensities[name] = intensity
        total_intensity += intensity
        equation_parts.append(f"10^({db:.2f}/10)")

    equation_str = " + ".join(equation_parts)
    
    # Calculate final total dB
    total_db = 10 * math.log10(total_intensity)

    print(f"\nFinal Equation:")
    print(f"Total dB = 10 * log10({equation_str})")
    print(f"Total dB = 10 * log10({total_intensity:.2f})")
    print(f"The total sound level you hear is: {total_db:.2f} dB")
    
    return total_db

# Execute the function
final_answer = calculate_total_sound_level()
# The final answer is wrapped according to the instruction.
# The print statements above provide the detailed breakdown.
# The final numerical answer is returned here.
# print(f"\n<<<{final_answer:.2f}>>>") # This line is for final answer extraction.

if __name__ == '__main__':
    # This block is for direct execution and will not be part of the response.
    # It helps in verifying the final answer.
    # The final answer is printed within the function itself.
    pass