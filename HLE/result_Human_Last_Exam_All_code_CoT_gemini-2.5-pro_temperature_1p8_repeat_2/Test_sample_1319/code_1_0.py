def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, besides Mount Olympus.
    """
    # Elevations are in meters.
    # Data is based on mountains explicitly mentioned in the text of the Iliad.
    mount_olympus_elevation = 2918

    # Dictionary of other historical mountains mentioned in the Iliad and their elevations
    iliad_mountains = {
        "Mount Ida (Turkey, overlooking Troy)": 1774,
        "Mount Samothrace (Poseidon's viewpoint)": 1611,
        "Mount Athos (on Hera's flight path)": 2033,
        "Mount Pelion (source of Achilles' spear)": 1624,
        "Mount Tmolus (home of the Maeonians)": 2159
    }

    # Find the tallest mountain from the dictionary
    tallest_mountain_name = ""
    max_elevation = 0
    for mountain, elevation in iliad_mountains.items():
        if elevation > max_elevation:
            max_elevation = elevation
            tallest_mountain_name = mountain
    
    print("Finding the tallest historical mountain in the Iliad (after Mount Olympus)...")
    print("-" * 50)
    print(f"Elevation of Mount Olympus: {mount_olympus_elevation} meters")
    print("Comparing other mountains mentioned in the epic:")
    for mountain, elevation in iliad_mountains.items():
        print(f"- {mountain}: {elevation} meters")
    
    print("-" * 50)
    print(f"The tallest mountain among these is {tallest_mountain_name}.")
    print(f"Final Answer: The tallest historical mountain mentioned in the Iliad after Mount Olympus is {tallest_mountain_name} with an elevation of {max_elevation} meters.")

if __name__ == "__main__":
    find_tallest_mountain_in_iliad()
