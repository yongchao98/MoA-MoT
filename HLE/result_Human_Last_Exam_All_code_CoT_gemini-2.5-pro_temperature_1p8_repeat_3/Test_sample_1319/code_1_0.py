def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # Dictionary of mountains mentioned in the Iliad (excluding Olympus) and their elevations in meters.
    mountains = {
        "Mount Tmolus": 2157,
        "Mount Athos": 2033,
        "Mount Ida (Trojan)": 1774,
        "Mount Othrys": 1726,
        "Mount Pelion": 1624,
        "Samothrace (Mount Fengari)": 1611,
    }

    # Find the tallest mountain from the dictionary
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print("Comparing the heights of historical mountains mentioned in the Iliad (after Mount Olympus):")
    
    # Print each mountain and its height for comparison
    # This loop demonstrates the "final equation" by showing all the numbers being compared.
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    print(f"\nThe tallest among these is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

if __name__ == "__main__":
    find_tallest_mountain_in_iliad()