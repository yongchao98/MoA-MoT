def find_tallest_mountain():
    """
    Finds the tallest historical mountain mentioned in the Iliad, besides Mount Olympus.
    """
    # A dictionary of historical mountains from the Iliad and their elevations in meters.
    # Mount Olympus is excluded as per the prompt.
    mountains = {
        "Mount Ida (Kaz Dağı)": 1774,
        "Mount Samothrace (Fengari)": 1611,
        "Mount Athos": 2033,
        "Mount Pelion": 1624,
        "Mount Ossa": 1978,
    }

    print("Comparing the elevations of mountains mentioned in the Iliad (other than Olympus):\n")

    # Initialize variables to store the tallest mountain found so far
    tallest_mountain_name = ""
    max_elevation = 0

    # Print each mountain and its elevation, and find the maximum
    for name, elevation in mountains.items():
        print(f"{name}: {elevation} meters")
        if elevation > max_elevation:
            max_elevation = elevation
            tallest_mountain_name = name

    print(f"\nAfter comparing the heights, the tallest mountain is {tallest_mountain_name} with an elevation of {max_elevation} meters.")

find_tallest_mountain()