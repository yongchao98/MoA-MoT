import operator

def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, and prints the result and comparison.
    """
    # Mountains mentioned in the Iliad (excluding Olympus) and their elevations in meters
    mountains = {
        "Mount Ida (Turkey)": 1774,
        "Mount Samothrace (Greece)": 1611,
        "Mount Athos (Greece)": 2033,
        "Mount Pelion (Greece)": 1624
    }

    # Find the tallest mountain
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name} at {tallest_mountain_height} meters.")
    print("\nHere is a comparison of the mountain elevations:")

    # Sort the mountains by height in descending order for the comparison "equation"
    sorted_mountains = sorted(mountains.items(), key=operator.itemgetter(1), reverse=True)
    
    # Build and print the comparison string, showing each name and number
    comparison_equation = []
    for name, height in sorted_mountains:
        comparison_equation.append(f"{name} ({height}m)")

    print(" > ".join(comparison_equation))

if __name__ == "__main__":
    find_tallest_iliad_mountain()