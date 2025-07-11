def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus is excluded as per the request.
    mountains = {
        "Mount Ida": 1774,
        "Mount Samothrace": 1611,
        "Mount Athos": 2033,
        "Mount Pelion": 1624
    }

    # Find the mountain with the maximum elevation
    tallest_mountain_name = ""
    max_height = 0
    for mountain, height in mountains.items():
        if height > max_height:
            max_height = height
            tallest_mountain_name = mountain

    # Print the comparison
    print("Comparing the heights of mountains from the Iliad (excluding Olympus):")
    print(f"Mount Ida: {mountains['Mount Ida']} m")
    print(f"Mount Samothrace: {mountains['Mount Samothrace']} m")
    print(f"Mount Athos: {mountains['Mount Athos']} m")
    print(f"Mount Pelion: {mountains['Mount Pelion']} m")
    print("\n----------------------------------------------------")
    # Print the final result
    print(f"The tallest among these historical mountains is {tallest_mountain_name} at {max_height} meters.")


if __name__ == "__main__":
    find_tallest_iliad_mountain()