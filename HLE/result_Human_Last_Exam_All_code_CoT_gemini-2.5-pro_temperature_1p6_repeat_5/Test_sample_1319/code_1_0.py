def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, besides Mount Olympus.
    """
    # Dictionary of mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus (2918m) is excluded as per the prompt.
    iliad_mountains = {
        "Mount Ida (near Troy)": 1774,
        "Mount Samothrace": 1611,
        "Mount Pelion": 1624,
        "Mount Athos": 2033
    }

    # Find the mountain with the maximum elevation
    tallest_mountain_name = max(iliad_mountains, key=iliad_mountains.get)
    tallest_mountain_height = iliad_mountains[tallest_mountain_name]

    print(f"The mountains mentioned in the Iliad (besides Olympus) and their elevations are:")
    for name, height in iliad_mountains.items():
        print(f"- {name}: {height} meters")

    print(f"\nThe tallest among these is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

find_tallest_iliad_mountain()