def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of historical mountains from the Iliad and their elevations in meters.
    # Mount Olympus (2918m) is excluded as per the prompt.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida (near Troy)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace (Fengari)": 1611
    }

    print("Comparing the elevations of historical mountains mentioned in the Iliad (excluding Mount Olympus):\n")
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    # Find the tallest mountain from the dictionary
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    # Create the comparison equation string with all numbers
    # Sort mountains by height in descending order to build the equation
    sorted_heights = sorted(mountains.values(), reverse=True)
    
    # Format the numbers into a comparison string "2033 > 1774 > ..."
    equation_str = " > ".join(map(str, sorted_heights))

    print(f"\nThe comparison of heights is: {equation_str}")

    print(f"\nTherefore, after Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")


if __name__ == "__main__":
    find_tallest_iliad_mountain()