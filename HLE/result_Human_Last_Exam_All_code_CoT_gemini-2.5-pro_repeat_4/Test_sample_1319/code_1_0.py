def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    after Mount Olympus, by comparing their known elevations.
    """
    # A dictionary of mountains mentioned in the Iliad (excluding Olympus)
    # and their respective elevations in meters.
    mountains = {
        "Mount Ida (Troad)": 1774,
        "Mount Samothrace": 1611,
        "Mount Athos": 2033,
        "Mount Pelion": 1624
    }

    # Find the tallest mountain from the dictionary
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_elevation = mountains[tallest_mountain_name]

    print("Comparing the elevations of historical mountains from the Iliad (excluding Olympus):")
    for name, elevation in mountains.items():
        print(f"- {name}: {elevation} meters")

    print("\nThe comparison can be represented by finding the maximum value:")
    # The prompt asks to "output each number in the final equation!"
    # Here we show the numbers being compared to find the maximum.
    all_elevations = list(mountains.values())
    print(f"max({all_elevations[0]}, {all_elevations[1]}, {all_elevations[2]}, {all_elevations[3]}) = {tallest_mountain_elevation}")

    print(f"\nThe tallest mountain among these is {tallest_mountain_name}.")


if __name__ == "__main__":
    find_tallest_iliad_mountain()
<<<Mount Athos>>>