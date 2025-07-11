def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of historical mountains mentioned in the Iliad (excluding Olympus)
    # and their respective modern elevations in meters.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    # Find the mountain with the highest elevation
    tallest_mountain_name = ""
    max_elevation = 0
    for name, elevation in mountains.items():
        if elevation > max_elevation:
            max_elevation = elevation
            tallest_mountain_name = name

    # Print the comparison
    print("Comparing the elevations of mountains mentioned in the Iliad (besides Olympus):")
    for name, elevation in mountains.items():
        print(f"- {name}: {elevation}m")

    # Print the equation as requested
    elevations = list(mountains.values())
    print(f"\nThe final equation to find the maximum height is: max({elevations[0]}, {elevations[1]}, {elevations[2]}, {elevations[3]}) = {max_elevation}")

    # Print the final answer
    print(f"\nAfter Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name} with an elevation of {max_elevation} meters.")

if __name__ == "__main__":
    find_tallest_iliad_mountain()