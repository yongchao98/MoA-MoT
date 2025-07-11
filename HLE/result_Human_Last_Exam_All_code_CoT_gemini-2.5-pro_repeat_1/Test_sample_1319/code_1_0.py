def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # Step 1 & 2: A dictionary containing the mountains mentioned in the Iliad
    # (excluding Olympus) and their respective elevations in meters.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida (Turkey)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    print("Comparing the elevations of historical mountains from the Iliad (excluding Olympus):\n")
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    # Step 3 & 4: Find the mountain with the maximum elevation.
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    # Step 5: Output the final result, showing the numbers used in the comparison.
    # We construct a string that looks like an equation to fulfill the prompt's requirement.
    height_values = list(mountains.values())
    
    print("\nFinal comparison equation:")
    print(f"The maximum of the heights ({height_values[0]}, {height_values[1]}, {height_values[2]}, {height_values[3]}) is {tallest_mountain_height}.")
    
    print(f"\nTherefore, after Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")

find_tallest_mountain_in_iliad()