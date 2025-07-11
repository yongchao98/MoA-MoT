import math

def solve_sound_level():
    """
    Calculates the total sound level at a new position from multiple sources.
    """
    # 1. Establish coordinates and new position
    # My initial position is (0,0)
    # The construction is in front (positive y), so its direction is (0, 1).
    # I walk 25m towards it.
    my_new_pos = (0, 25)

    # 2. Define sound sources: (Name, Sound Pressure Level at 1m, position (x, y))
    sources = [
        ("Dog", 55, (-25, 0)),        # 25m to the left
        ("Train", 110, (50, 0)),        # 50m to the right
        ("Construction", 90, (0, 75)),        # 75m in front
        ("People", 75, (0, -10))      # 10m behind
    ]

    print(f"My new position is {my_new_pos}\n")
    print("--- Calculating sound level from each source at the new position ---\n")

    # A list to store the calculated sound levels at my new position
    sound_levels_at_my_pos = []

    # 3. Calculate new distances and sound levels for each source
    for name, spl_at_1m, pos in sources:
        # Calculate the distance from my new position to the source
        distance = math.sqrt((pos[0] - my_new_pos[0])**2 + (pos[1] - my_new_pos[1])**2)

        # Calculate the sound level (Lp2) at this new distance
        # Lp2 = Lp1 - 20 * log10(r2)
        spl_at_new_pos = spl_at_1m - 20 * math.log10(distance)
        sound_levels_at_my_pos.append(spl_at_new_pos)

        print(f"{name}:")
        print(f"  - New distance: {distance:.2f} meters")
        print(f"  - Sound level: {spl_at_1m} - 20 * log10({distance:.2f}) = {spl_at_new_pos:.2f} dB\n")

    # 4. Combine the sound levels
    # Convert each dB to a value proportional to intensity (10^(Lp/10)), sum them,
    # and convert back to dB.
    total_intensity_proportional = sum([10**(spl / 10) for spl in sound_levels_at_my_pos])
    total_spl = 10 * math.log10(total_intensity_proportional)

    # 5. Print the final equation and result
    print("--- Calculating the total combined sound level ---")
    print("To combine the sounds, we convert each decibel value to an intensity, sum them, and convert back to decibels.\n")

    equation_parts = [f"10^({spl:.2f}/10)" for spl in sound_levels_at_my_pos]
    final_equation = "10 * log10( " + " + ".join(equation_parts) + " )"

    print("Final Calculation:")
    print(f"Total Sound Level = {final_equation}\n")
    print(f"The total sound level I hear is {total_spl:.2f} dB.")

solve_sound_level()
<<<75.14>>>