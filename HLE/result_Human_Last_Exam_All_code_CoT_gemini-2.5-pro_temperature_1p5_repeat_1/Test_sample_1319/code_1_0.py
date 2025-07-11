def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    The data is pre-researched and stored in a dictionary.
    """
    mountains = {
        "Mount Ida (Crete)": 2456,
        "Mount Athos": 2033,
        "Mount Ossa": 1978,
        "Mount Ida (Turkey)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    # Find the mountain with the maximum elevation
    tallest_mountain_name = ""
    max_elevation = 0
    for name, elevation in mountains.items():
        if elevation > max_elevation:
            max_elevation = elevation
            tallest_mountain_name = name

    print("Comparing the elevations of historical mountains from the Iliad (excluding Olympus):")
    # Print each mountain and its elevation
    # The 'equation' here is the comparison of all listed numbers to find the maximum.
    for name, elevation in mountains.items():
        print(f"- {name}: {elevation} meters")

    print("\nAfter comparing the values:")
    print(f"2456 (Mount Ida, Crete) > 2033 (Mount Athos) > 1978 (Mount Ossa) > 1774 (Mount Ida, Turkey) > 1624 (Mount Pelion) > 1611 (Mount Samothrace)")
    print(f"\nThe tallest mountain is {tallest_mountain_name} with an elevation of {max_elevation} meters.")

if __name__ == "__main__":
    find_tallest_mountain()