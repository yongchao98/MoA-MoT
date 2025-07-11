def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, by comparing their elevations.
    """
    # Dictionary of historical mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus (2917 m) is excluded as per the user's request.
    mountains = {
        "Mount Ida (near Troy)": 1774,
        "Mount Samothrace": 1611,
        "Mount Athos": 2033,
        "Mount Pelion": 1624
    }

    print("Comparing the elevations of historical mountains mentioned in the Iliad (besides Olympus):\n")

    tallest_mountain_name = ""
    max_height = 0

    # Iterate through the dictionary to find the tallest mountain
    # and print each mountain's data as part of the process.
    for name, height in mountains.items():
        # This print statement fulfills the requirement to "output each number".
        print(f"- {name}: {height} meters")
        if height > max_height:
            max_height = height
            tallest_mountain_name = name

    print(f"\nAfter comparing the elevations, the tallest historical mountain mentioned is {tallest_mountain_name} with a height of {max_height} meters.")

find_tallest_mountain()