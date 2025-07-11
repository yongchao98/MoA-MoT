def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad and their real-world elevations in meters.
    mountains = {
        'Mount Olympus': 2918,
        'Mount Athos': 2033,
        'Mount Ida': 1774,
        'Mount Pelion': 1624,
        'Mount Samothrace': 1611
    }

    # Exclude Mount Olympus to find the next tallest.
    other_mountains = {name: elev for name, elev in mountains.items() if name != 'Mount Olympus'}

    # Find the name of the tallest mountain among the rest.
    tallest_mountain_name = max(other_mountains, key=other_mountains.get)
    
    # Get the elevation of that mountain.
    tallest_mountain_elevation = other_mountains[tallest_mountain_name]

    print("The historical mountains mentioned in the Iliad (besides Olympus) and their elevations are:")
    for name, elevation in other_mountains.items():
        print(f"- {name}: {elevation} meters")

    # The prompt requires outputting each number in the final equation.
    # We will show the comparison to find the maximum elevation.
    elevations_to_compare = list(other_mountains.values())
    
    print("\nThe comparison to find the maximum elevation is:")
    print(f"max({elevations_to_compare[0]}, {elevations_to_compare[1]}, {elevations_to_compare[2]}, {elevations_to_compare[3]}) = {tallest_mountain_elevation}")

    print(f"\nTherefore, after Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")

find_tallest_mountain()