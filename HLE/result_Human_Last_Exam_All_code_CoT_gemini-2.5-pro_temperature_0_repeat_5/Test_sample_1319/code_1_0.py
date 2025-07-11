def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of historical mountains mentioned by name in the Iliad
    # and their elevations in meters. Mount Olympus (2918m) is excluded as per the prompt.
    iliad_mountains = {
        "Mount Tmolus": 2157,
        "Mount Athos": 2033,
        "Mount Ida (Trojan)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    print("The historical mountains mentioned by name in the Iliad (excluding Mount Olympus) are:")
    # Print each mountain and its elevation, which satisfies the "output each number" requirement.
    for mountain, elevation in iliad_mountains.items():
        print(f"- {mountain}: {elevation} meters")

    # Find the mountain with the highest elevation in the dictionary
    if not iliad_mountains:
        print("\nNo mountains to compare.")
        return

    tallest_mountain_name = max(iliad_mountains, key=iliad_mountains.get)
    tallest_mountain_elevation = iliad_mountains[tallest_mountain_name]

    print(f"\nBy comparing the elevations, the tallest mountain in this list is {tallest_mountain_name}.")
    print(f"\nTherefore, after Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name} with an elevation of {tallest_mountain_elevation} meters.")

if __name__ == "__main__":
    find_tallest_mountain_in_iliad()