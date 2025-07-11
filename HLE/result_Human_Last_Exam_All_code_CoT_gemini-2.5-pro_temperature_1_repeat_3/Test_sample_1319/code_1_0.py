import operator

def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus,
    and prints the result along with a comparison.
    """
    # A dictionary of historical mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus (2918m) is excluded as per the user's request.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida (near Troy)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Sort the mountains based on their elevation in descending order.
    # The key for sorting is the second item (index 1) of each dictionary item, which is the elevation.
    sorted_mountains = sorted(mountains.items(), key=operator.itemgetter(1), reverse=True)

    # The tallest mountain is the first one in the sorted list.
    tallest_mountain_name = sorted_mountains[0][0]
    tallest_mountain_elevation = sorted_mountains[0][1]

    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")

    # Create the comparison string, which serves as the "equation" showing each mountain and its height.
    comparison_parts = []
    for name, elevation in sorted_mountains:
        # Each part is a string like "Mount Athos (2033m)"
        comparison_parts.append(f"{name} ({elevation}m)")

    # Join the parts with " > " to form the final comparison equation.
    comparison_equation = " > ".join(comparison_parts)

    print("\nHere is the height comparison:")
    print(comparison_equation)

if __name__ == "__main__":
    find_tallest_mountain_in_iliad()