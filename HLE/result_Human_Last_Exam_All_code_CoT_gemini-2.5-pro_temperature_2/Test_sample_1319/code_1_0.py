def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    The elevations are the modern accepted heights in meters.
    """
    
    # Dictionary of historical mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus (2918m) is excluded as per the request.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida (Kaz Dağı)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace (Fengari)": 1611
    }

    print("Comparing the heights of mountains mentioned in the Iliad (excluding Olympus):\n")
    
    # Print each mountain and its height. This fulfills the "output each number" requirement.
    for name, height in mountains.items():
        print(f"{name}: {height} meters")

    # Find the mountain with the maximum height
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print(f"\nAfter comparing the elevations, the tallest historical mountain mentioned in the Iliad, after Mount Olympus, is {tallest_mountain_name}.")

if __name__ == "__main__":
    find_tallest_iliad_mountain()
