def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, and prints the result.
    """
    # A dictionary of mountains mentioned in the Iliad and their approximate elevations in meters.
    mountains_in_iliad = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Create a new dictionary excluding Mount Olympus.
    other_mountains = {name: height for name, height in mountains_in_iliad.items() if name != "Mount Olympus"}

    # Check if there are any other mountains to compare.
    if not other_mountains:
        print("No other mountains found to compare.")
        return

    # Find the name of the tallest mountain among the rest.
    # The 'key' argument for max() specifies a function to be called on each dictionary key
    # to get the value to compare. Here, other_mountains.get returns the elevation.
    tallest_mountain_name = max(other_mountains, key=other_mountains.get)
    tallest_mountain_height = other_mountains[tallest_mountain_name]

    print("After Mount Olympus, the historical mountains in the Iliad are ranked by height as follows:")
    # Sort the remaining mountains by height in descending order for printing.
    sorted_mountains = sorted(other_mountains.items(), key=lambda item: item[1], reverse=True)

    for name, height in sorted_mountains:
        print(f"- {name} has an elevation of {height} meters.")

    print(f"\nThe tallest historical mountain mentioned in the Iliad after Mount Olympus is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")


if __name__ == "__main__":
    find_tallest_mountain()