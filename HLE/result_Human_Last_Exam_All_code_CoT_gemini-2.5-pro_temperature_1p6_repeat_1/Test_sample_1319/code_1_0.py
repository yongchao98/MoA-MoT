import operator

def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    after Mount Olympus, by comparing their elevations.
    """
    # A dictionary of mountains mentioned in the Iliad (excluding Olympus)
    # and their respective elevations in meters.
    mountains = {
        "Mount Ida (Turkey)": 1774,
        "Mount Samothrace (Greece)": 1611,
        "Mount Athos (Greece)": 2033,
        "Mount Pelion (Greece)": 1624,
        "Mount Ossa (Greece)": 1978,
    }

    print("Comparing the elevations of mountains mentioned in the Iliad (after Olympus):")
    
    # Print each mountain and its elevation
    # This fulfills the request to show the numbers used in the "equation" or comparison.
    for mountain, elevation in mountains.items():
        print(f"- {mountain}: {elevation} meters")

    # Find the mountain with the highest elevation
    if not mountains:
        print("\nNo mountains to compare.")
        return

    # Using max with a key to find the entry with the highest value
    tallest_mountain_name, tallest_mountain_elevation = max(mountains.items(), key=operator.itemgetter(1))

    print(f"\nThe tallest among these is {tallest_mountain_name} with an elevation of {tallest_mountain_elevation} meters.")


if __name__ == "__main__":
    find_tallest_iliad_mountain()
