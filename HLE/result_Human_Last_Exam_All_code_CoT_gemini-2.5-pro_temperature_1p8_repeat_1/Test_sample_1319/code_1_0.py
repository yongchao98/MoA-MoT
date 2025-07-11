def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # Step 1 & 2: Store mountains and their elevations (in meters).
    mountains = {
        'Mount Athos': 2033,
        'Mount Ida': 1774,
        'Mount Pelion': 1624,
        'Mount Samothrace': 1611
    }

    # Step 3: Find the mountain with the maximum elevation.
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_elevation = mountains[tallest_mountain_name]

    # Step 4: Output the comparison and the result.
    # The prompt requires outputting each number in the final equation.
    print("Comparing the elevations of mountains from the Iliad (excluding Olympus):")
    
    # Get all elevation values for the "equation" part.
    elevations = list(mountains.values())
    
    print(f"The elevations considered are: {mountains['Mount Athos']}m (Athos), {mountains['Mount Ida']}m (Ida), {mountains['Mount Pelion']}m (Pelion), and {mountains['Mount Samothrace']}m (Samothrace).")
    
    # Printing the "equation" as a comparison of all numbers
    print(f"\nThe comparison equation is: max({elevations[0]}, {elevations[1]}, {elevations[2]}, {elevations[3]}) = {tallest_mountain_elevation}")
    
    print(f"\nAfter Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")

# Execute the function
find_tallest_iliad_mountain()
