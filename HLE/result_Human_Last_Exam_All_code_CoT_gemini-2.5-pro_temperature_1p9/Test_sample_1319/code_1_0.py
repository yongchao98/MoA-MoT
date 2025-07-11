def find_tallest_mountain():
    """
    This script identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, by comparing the elevations of other mentioned mountains.
    """
    # A dictionary of mountains mentioned in the Iliad (besides Olympus) and their elevations in meters.
    # The elevations are based on their modern-day identities.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida (Kaz Dağı, Turkey)": 1774,
        "Samothrace (Mount Fengari)": 1611,
        "Mount Pelion": 1624,
    }

    print("Comparing the elevations of historical mountains mentioned in the Iliad (excluding Olympus):\n")
    for mountain, elevation in mountains.items():
        print(f"- {mountain}: {elevation} meters")

    # Find the mountain with the highest elevation
    tallest_mountain_name = ""
    tallest_mountain_elevation = 0
    for mountain, elevation in mountains.items():
        if elevation > tallest_mountain_elevation:
            tallest_mountain_elevation = elevation
            tallest_mountain_name = mountain

    print(f"\nThe comparison shows that the tallest among these is {tallest_mountain_name}.")
    print(f"\nFinal Equation: max(2033, 1774, 1611, 1624) = 2033")
    print(f"\nTherefore, after Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name} with an elevation of {tallest_mountain_elevation} meters.")

find_tallest_mountain()