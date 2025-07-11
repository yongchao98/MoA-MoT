import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1: Define initial sound levels at 1 meter and source positions.
    # The format is {name: {"L1": dB at 1m, "pos": (x, y)}}.
    sources = {
        "Dog": {"L1": 55, "pos": (-25, 0)},
        "Train": {"L1": 110, "pos": (50, 0)},
        "Construction": {"L1": 90, "pos": (0, 75)},
        "People": {"L1": 75, "pos": (0, -10)}
    }
    r1 = 1.0  # Reference distance is 1 meter.

    # Step 2: Define the observer's new position.
    my_pos = (0, 25)

    total_intensity = 0
    new_levels_db = {}

    print("Calculating individual sound levels at the new position (0, 25):\n")

    # Step 3 & 4: Calculate new distance and new sound level for each source.
    for name, data in sources.items():
        source_pos = data["pos"]
        L1 = data["L1"]

        # Calculate the new distance (r2) from the observer to the source.
        r2 = math.sqrt((source_pos[0] - my_pos[0])**2 + (source_pos[1] - my_pos[1])**2)

        # Calculate the new sound level (L2) at distance r2.
        # Formula: L2 = L1 - 20 * log10(r2 / r1). Since r1=1, this is L1 - 20 * log10(r2).
        L2 = L1 - 20 * math.log10(r2)
        new_levels_db[name] = L2

        print(f"Distance to {name}: {r2:.2f} m, New Sound Level: {L2:.2f} dB")

        # Step 5 (part 1): Convert the new sound level (L2) to intensity and add to the total.
        # Formula: I = 10^(L/10) (using relative intensity).
        intensity = 10**(L2 / 10)
        total_intensity += intensity

    # Step 5 (part 2) & 6: Convert total intensity back to a total sound level.
    # Formula: L_total = 10 * log10(I_total)
    L_total = 10 * math.log10(total_intensity)

    print("\n--------------------------------------------------------------")
    print("The final combined sound level is calculated from the sum of the intensities.")
    print("Equation: L_total = 10 * log10( 10^(L_dog/10) + 10^(L_train/10) + 10^(L_construction/10) + 10^(L_people/10) )")
    print("Plugging in the numbers for the new sound levels:")
    
    # Show the final equation with the calculated values.
    dog_level = new_levels_db["Dog"]
    train_level = new_levels_db["Train"]
    construction_level = new_levels_db["Construction"]
    people_level = new_levels_db["People"]

    print(f"L_total = 10 * log10( 10^({dog_level:.2f}/10) + 10^({train_level:.2f}/10) + 10^({construction_level:.2f}/10) + 10^({people_level:.2f}/10) )")
    print("--------------------------------------------------------------\n")
    print(f"The total sound level at your new location is: {L_total:.2f} dB")


calculate_total_sound_level()
<<<75.14>>>