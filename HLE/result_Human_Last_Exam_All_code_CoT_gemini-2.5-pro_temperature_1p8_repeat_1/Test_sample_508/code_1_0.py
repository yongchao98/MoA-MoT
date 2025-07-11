import math

def calculate_total_sound_level():
    """
    Calculates the total sound level at a specific location from multiple sources.
    """
    # 1. Define source data: (name, (x, y), dB @ 1m)
    sources = [
        ("Dog", (-25, 0), 55),
        ("Train", (50, 0), 110),
        ("Construction", (0, 75), 90),
        ("People", (0, -10), 75)
    ]

    # 2. My new position: walked 25m towards the construction site (0, 75) from (0,0)
    my_new_pos = (0, 25)

    print(f"Your new location is at coordinates {my_new_pos}.")
    print("-" * 30)

    # 3. Calculate new sound levels for each source
    new_db_levels = []
    sum_of_intensities_ratio = 0
    equation_parts = []

    for name, source_pos, db_at_1m in sources:
        # Calculate new distance
        distance = math.sqrt((source_pos[0] - my_new_pos[0])**2 + (source_pos[1] - my_new_pos[1])**2)
        
        # Calculate new sound level at the new distance
        # Formula: L2 = L1 - 20 * log10(r2), since r1=1m
        new_db = db_at_1m - 20 * math.log10(distance)
        new_db_levels.append(new_db)
        
        print(f"Source: {name}")
        print(f"New distance: {distance:.2f} meters")
        print(f"Sound level from this source: {new_db:.2f} dB")
        print("-" * 30)

        # Prepare for total calculation
        sum_of_intensities_ratio += 10**(new_db / 10)
        equation_parts.append(f"10^({new_db:.2f}/10)")
        
    # 4. Calculate total sound level by summing the intensity ratios
    total_db = 10 * math.log10(sum_of_intensities_ratio)

    # 5. Print the final result and equation
    print("\nTo find the total sound level, we use the formula:")
    print("Total dB = 10 * log10( Σ(10^(L/10)) )\n")
    
    final_equation_str = " + ".join(equation_parts)
    print("Calculation:")
    print(f"Total Sound Level = 10 * log10({final_equation_str})")
    print(f"Total Sound Level = {total_db:.2f} dB")


if __name__ == "__main__":
    calculate_total_sound_level()
    # The final answer needs to be extracted for the <<<...>>> format.
    # Performing the calculation again here just to store the final value.
    sources = [((-25, 0), 55), ((50, 0), 110), ((0, 75), 90), ((0, -10), 75)]
    my_new_pos = (0, 25)
    sum_of_intensities_ratio = 0
    for source_pos, db_at_1m in sources:
        distance = math.sqrt((source_pos[0] - my_new_pos[0])**2 + (source_pos[1] - my_new_pos[1])**2)
        new_db = db_at_1m - 20 * math.log10(distance)
        sum_of_intensities_ratio += 10**(new_db / 10)
    total_db = 10 * math.log10(sum_of_intensities_ratio)
    # The required output format
    print(f"\n<<<{total_db:.2f}>>>")
