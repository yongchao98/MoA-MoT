def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    The script stores mountain data, performs a comparison, and prints the result.
    """
    # Dictionary of mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus is included for context but will be excluded from the main comparison.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Samothrace": 1611,
    }

    # Exclude Mount Olympus for the comparison as requested.
    contender_mountains = {name: height for name, height in mountains.items() if name != "Mount Olympus"}

    # Find the tallest mountain among the contenders.
    if not contender_mountains:
        print("No mountains to compare after excluding Mount Olympus.")
        return

    tallest_mountain_name = max(contender_mountains, key=contender_mountains.get)
    tallest_mountain_height = contender_mountains[tallest_mountain_name]

    # Create a simple "equation" like comparison string showing all numbers.
    sorted_heights = sorted(contender_mountains.values(), reverse=True)
    comparison_equation = " > ".join(map(str, sorted_heights))

    print("After Mount Olympus, the notable historical mountains in the Iliad are compared.")
    print(f"The comparison of their heights is: {tallest_mountain_name} ({tallest_mountain_height}m) is the tallest.")
    print(f"The equation representing the comparison of elevations is: {comparison_equation}")
    print(f"\nTherefore, the tallest historical mountain mentioned in the Iliad after Mount Olympus is {tallest_mountain_name}.")


if __name__ == '__main__':
    find_tallest_mountain_in_iliad()