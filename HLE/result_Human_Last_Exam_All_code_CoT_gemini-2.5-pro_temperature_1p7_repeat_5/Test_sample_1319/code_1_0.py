def find_tallest_mountain():
    """
    This function finds the tallest historical mountain mentioned in the Iliad,
    other than Mount Olympus.
    """
    # A dictionary of mountains from the Iliad and their elevations in meters.
    # Mount Olympus (2917m) is excluded as per the request.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Find the tallest mountain by finding the maximum value in the dictionary.
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_elevation = mountains[tallest_mountain_name]

    # Print the final answer and the comparison.
    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")
    print(f"Its elevation is {tallest_mountain_elevation} meters.")

    print("\nElevation Comparison:")
    # Sort mountains by elevation in descending order for clear comparison.
    for name, elevation in sorted(mountains.items(), key=lambda item: item[1], reverse=True):
        print(f"- {name}: {elevation} m")


if __name__ == "__main__":
    find_tallest_mountain()