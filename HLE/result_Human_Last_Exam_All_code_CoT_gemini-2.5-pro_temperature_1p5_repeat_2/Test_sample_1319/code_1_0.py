def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, besides Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad (excluding Olympus) and their elevations in meters.
    # Sources for mountains: Homer's Iliad.
    # Sources for elevations: Modern geographical data.
    mountains = {
        "Mount Ida": 1774,      # Where Zeus and other gods watched the Trojan War.
        "Mount Samothrace": 1611, # Where Poseidon watched the war.
        "Mount Athos": 2033,      # Mentioned in Hera's journey.
        "Mount Pelion": 1624      # The source of the ash wood for Achilles' spear.
    }

    # Find the tallest mountain from the dictionary
    tallest_mountain_name = ""
    tallest_mountain_height = 0
    for name, height in mountains.items():
        if height > tallest_mountain_height:
            tallest_mountain_height = height
            tallest_mountain_name = name

    print("Comparing the heights of mountains from the Iliad (excluding Olympus):")
    # Using a formatted print to show each number in the final comparison
    print(f" - Mount Athos: {mountains['Mount Athos']} meters")
    print(f" - Mount Pelion: {mountains['Mount Pelion']} meters")
    print(f" - Mount Ida: {mountains['Mount Ida']} meters")
    print(f" - Mount Samothrace: {mountains['Mount Samothrace']} meters")
    print("\nAfter comparing the elevations, the tallest historical mountain mentioned in the Iliad besides Mount Olympus is:")
    print(f"{tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")


if __name__ == '__main__':
    find_tallest_mountain_in_iliad()