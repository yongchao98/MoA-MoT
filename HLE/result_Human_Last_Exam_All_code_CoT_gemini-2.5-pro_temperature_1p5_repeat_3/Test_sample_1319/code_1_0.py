def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # Dictionary of mountains from the Iliad (excluding Olympus) and their elevations in meters.
    # These are the actual historical mountains.
    mountains = {
        "Mount Ida": 1774,
        "Mount Samothrace": 1611,
        "Mount Athos": 2033,
        "Mount Pelion": 1624
    }

    # Find the mountain with the maximum elevation
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print("Comparing the elevations of historical mountains mentioned in the Iliad (besides Olympus):")
    
    # "Equation" part: Printing each mountain and its height
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    print("\nBased on the data, the tallest mountain is:")
    print(f"{tallest_mountain_name} at {tallest_mountain_height} meters.")

find_tallest_mountain()