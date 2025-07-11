def find_tallest_iliad_mountain():
    """
    Identifies the tallest mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # A dictionary of prominent mountains from the Iliad and their elevations in meters.
    mountains = {
        'Olympus': 2918,
        'Athos': 2033,
        'Ida': 1774,
        'Pelion': 1624,
        'Samothrace': 1611
    }

    # As per the request, we exclude Mount Olympus from our search.
    # We create a new dictionary without it.
    contenders = {name: height for name, height in mountains.items() if name != 'Olympus'}

    # Find the mountain with the maximum elevation among the contenders.
    tallest_mountain_name = max(contenders, key=contenders.get)
    tallest_mountain_height = contenders[tallest_mountain_name]

    print(f"The mountains considered are: {list(contenders.keys())}")
    print("Comparing their heights:")
    for name, height in contenders.items():
        print(f"- {name}: {height} meters")

    print("\nAfter Mount Olympus, the tallest historical mountain mentioned in the Iliad is:")
    # Final answer output, including the name (text) and the elevation (number).
    print(f"{tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

if __name__ == '__main__':
    find_tallest_iliad_mountain()