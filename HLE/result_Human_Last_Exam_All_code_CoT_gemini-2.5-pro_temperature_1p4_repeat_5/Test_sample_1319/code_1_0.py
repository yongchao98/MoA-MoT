def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad (excluding Olympus) and their elevations in meters.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    print("Comparing the elevations of mountains mentioned in the Iliad (after Olympus):")

    tallest_mountain_name = ""
    max_elevation = 0

    # Iterate through the dictionary to print each mountain and find the tallest
    for name, elevation in mountains.items():
        print(f"- {name}: {elevation} meters")
        if elevation > max_elevation:
            max_elevation = elevation
            tallest_mountain_name = name

    print(f"\nThe tallest among these is {tallest_mountain_name} with an elevation of {max_elevation} meters.")

if __name__ == "__main__":
    find_tallest_mountain()