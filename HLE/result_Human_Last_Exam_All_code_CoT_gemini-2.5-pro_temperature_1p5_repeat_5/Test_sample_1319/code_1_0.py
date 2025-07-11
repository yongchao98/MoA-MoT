def solve_mountain_height_puzzle():
    """
    This script identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, by comparing their modern elevations.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus is included for completeness but will be excluded from the main comparison.
    mountains = {
        "Mount Olympus": 2917,
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    # Exclude Mount Olympus as per the user's request.
    historical_mountains = {name: height for name, height in mountains.items() if name != "Mount Olympus"}

    # Sort the mountains by elevation in descending order to find the tallest and create the equation.
    sorted_mountains = sorted(historical_mountains.items(), key=lambda item: item[1], reverse=True)

    # Find the tallest mountain among the historical ones.
    if sorted_mountains:
        tallest_mountain_name, tallest_mountain_height = sorted_mountains[0]

        print("Comparing the heights of historical mountains mentioned in the Iliad (excluding Mount Olympus):")

        # Build and print the comparison equation string.
        # This will output each mountain name and its height.
        equation_parts = [f"{name} ({height}m)" for name, height in sorted_mountains]
        equation_string = " > ".join(equation_parts)
        print(equation_string)
        print("\nBased on this, the tallest historical mountain mentioned in the Iliad after Mount Olympus is:")
        print(f"{tallest_mountain_name} at {tallest_mountain_height} meters.")
    else:
        print("No other historical mountains found for comparison.")

solve_mountain_height_puzzle()
