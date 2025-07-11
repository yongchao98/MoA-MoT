def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus,
    by comparing the elevations of the most prominent mountains featured in the epic.
    """
    # A dictionary of mountains from the Iliad and their real-world elevations in meters.
    # Mount Olympus is included for context but excluded from the final comparison.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida (Trojan)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Create a new dictionary excluding Mount Olympus.
    other_mountains = {name: elev for name, elev in mountains.items() if name != "Mount Olympus"}

    # Sort the mountains by height in descending order to prepare for the output.
    sorted_mountains = sorted(other_mountains.items(), key=lambda item: item[1], reverse=True)

    # The tallest mountain is the first one in the sorted list.
    tallest_name, tallest_height = sorted_mountains[0]

    print("Comparing the elevations of historical mountains from the Iliad (excluding Mount Olympus):")

    # Build and print the final comparison "equation", showing each number.
    equation_parts = []
    for name, height in sorted_mountains:
        # Add each number (height) and name to the equation.
        equation_parts.append(str(height))

    print("The comparison is: " + " > ".join(equation_parts))

    print(f"\nBased on this comparison, the tallest mountain after Mount Olympus is {tallest_name}.")

# Execute the function to print the result.
find_tallest_mountain_in_iliad()