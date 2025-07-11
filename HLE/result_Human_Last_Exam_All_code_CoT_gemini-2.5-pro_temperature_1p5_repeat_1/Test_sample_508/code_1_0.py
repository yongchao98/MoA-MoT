import math

def solve_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # 1. & 2. Define initial positions and new position
    my_new_pos = (0, 25)
    sources = {
        "Dog":          {"lp_at_1m": 55, "pos": (-25, 0)},
        "Train":        {"lp_at_1m": 110, "pos": (50, 0)},
        "Construction": {"lp_at_1m": 90, "pos": (0, 75)},
        "People":       {"lp_at_1m": 75, "pos": (0, -10)},
    }

    print("Step-by-step calculation of the total sound level:\n")

    calculated_levels = []
    total_intensity_term = 0

    # 3. & 4. Calculate new distances and sound levels for each source
    print("Calculating individual sound levels at the new location (0, 25):\n")
    for name, data in sources.items():
        # Calculate new distance from (0, 25)
        distance = math.sqrt((my_new_pos[0] - data["pos"][0])**2 + (my_new_pos[1] - data["pos"][1])**2)
        
        # Calculate sound level (Lp2) at the new distance
        # Formula: Lp2 = Lp(at 1m) - 20 * log10(distance)
        lp2 = data["lp_at_1m"] - 20 * math.log10(distance)
        calculated_levels.append(lp2)
        
        # Calculate the intensity term 10^(Lp2/10) for the final summation
        total_intensity_term += 10**(lp2 / 10)
        
        print(f"Source: {name}")
        print(f"  - New distance: {distance:.2f} m")
        print(f"  - Sound level at new location: {lp2:.2f} dB\n")

    # 5. Combine the sound levels
    total_db = 10 * math.log10(total_intensity_term)

    # Print the final equation with the numbers
    equation_parts = [f"10^({level:.2f}/10)" for level in calculated_levels]
    equation = " + ".join(equation_parts)
    
    print("--------------------------------------------------")
    print("Final Calculation:")
    print("The total sound level is calculated by summing the intensities of each source.\n")
    print(f"Total dB = 10 * log10( {equation} )")
    print(f"Total dB = 10 * log10( {total_intensity_term:.2f} )")
    print(f"\nThe total sound level you hear is: {total_db:.2f} dB")
    print("--------------------------------------------------")
    
    # Return final answer for the system
    return total_db

# Run the calculation and store the final answer
final_answer = solve_sound_level()

# Print final answer in the required format
# print(f"\n<<<{final_answer:.2f}>>>") # This is for final answer output for the platform.
# For this case, let's keep the script runnable for the user as is.
# The user can see the final printed result.
