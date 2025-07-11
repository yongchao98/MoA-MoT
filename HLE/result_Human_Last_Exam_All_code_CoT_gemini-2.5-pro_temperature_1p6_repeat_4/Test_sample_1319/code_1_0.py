def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # Mountains mentioned in the Iliad (excluding Olympus) and their elevations in meters.
    mountains = {
        "Mount Ida": 1774,
        "Mount Samothrace": 1611,
        "Mount Athos": 2033
    }

    print("Comparing the elevations of historical mountains mentioned in the Iliad (after Olympus):\n")

    # Print each mountain and its elevation for comparison
    for mountain, elevation in mountains.items():
        print(f"- {mountain}: {elevation} meters")

    # Find the tallest mountain from the dictionary
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_elevation = mountains[tallest_mountain_name]

    print(f"\nBased on the comparison, the tallest mountain after Olympus is {tallest_mountain_name} with an elevation of {tallest_mountain_elevation} meters.")

if __name__ == "__main__":
    find_tallest_mountain_in_iliad()