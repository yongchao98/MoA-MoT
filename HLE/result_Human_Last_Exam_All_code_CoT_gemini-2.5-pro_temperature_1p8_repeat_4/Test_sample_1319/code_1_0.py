def find_tallest_mountain_in_iliad():
    """
    This function identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, by comparing the elevations of other mentioned mountains.
    """
    # A dictionary mapping the mountains mentioned in the Iliad (besides Olympus)
    # to their approximate modern elevations in meters.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida (Kaz Dağı)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace (Fengari Peak)": 1611
    }

    print("Comparing the elevations of historical mountains from the Iliad (after Olympus):")
    # Print each mountain and its height for the comparison
    # This fulfills the requirement to "output each number in the final equation".
    for name, height in mountains.items():
        print(f"{name}: {height} meters")

    # Find the name of the mountain with the maximum height
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print(f"\nThe tallest among these mountains is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

if __name__ == "__main__":
    find_tallest_mountain_in_iliad()