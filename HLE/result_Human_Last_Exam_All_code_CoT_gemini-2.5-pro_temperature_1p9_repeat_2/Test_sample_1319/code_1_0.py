def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # A dictionary of prominent mountains from the Iliad and their real-world elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida (Troad)": 1774,
        "Mount Samothrace": 1611
    }

    print("Finding the tallest historical mountain mentioned in the Iliad (after Mount Olympus).\n")
    print("The historical mountains and their elevations are:")
    # Create a new dictionary without Mount Olympus for the comparison
    historical_mountains_for_comparison = {name: height for name, height in mountains.items() if name != "Mount Olympus"}

    # Print each mountain being compared
    for name, height in historical_mountains_for_comparison.items():
        print(f"- {name}: {height} m")

    # Find the mountain with the maximum height from the filtered list
    tallest_mountain_name = max(historical_mountains_for_comparison, key=historical_mountains_for_comparison.get)
    tallest_mountain_height = historical_mountains_for_comparison[tallest_mountain_name]

    # Forming the "equation" or comparison string
    comparison_values = sorted(historical_mountains_for_comparison.values(), reverse=True)
    comparison_string = " > ".join(map(str, comparison_values))
    
    print(f"\nThe height comparison is: {comparison_string}")
    print(f"\nThe tallest mountain after Mount Olympus is {tallest_mountain_name}.")


if __name__ == "__main__":
    find_tallest_mountain()