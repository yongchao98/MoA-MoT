def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    after Mount Olympus, and prints the comparison.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    # Mt. Ida is the one near Troy (Kaz Dağı), not the one on Crete.
    mountains = {
        "Olympus": 2918,
        "Athos": 2033,
        "Ossa": 1978,
        "Ida": 1774,
        "Pelion": 1624,
        "Samothrace": 1611
    }

    # Exclude Mount Olympus from the comparison
    olympus_elevation = mountains.pop("Olympus")
    print(f"Excluding Mount Olympus, which has an elevation of {olympus_elevation} meters.\n")

    # Find the tallest mountain among the remaining ones
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_elevation = mountains[tallest_mountain_name]

    print("Comparing the elevations of the other mountains mentioned in the Iliad:")
    
    # Print each mountain and its elevation for the comparison "equation"
    for name, elevation in mountains.items():
        print(f"- {name}: {elevation} meters")

    print(f"\nThe tallest mountain after Olympus is {tallest_mountain_name} with an elevation of {tallest_mountain_elevation} meters.")

if __name__ == "__main__":
    find_tallest_mountain()