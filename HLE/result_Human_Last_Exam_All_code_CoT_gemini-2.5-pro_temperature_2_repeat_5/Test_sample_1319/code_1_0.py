def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    after Mount Olympus, and prints a comparison.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    mountains = {
        'Mount Olympus': 2918,
        'Mount Athos': 2033,
        'Mount Ida': 1774,
        'Mount Pelion': 1624,
        'Mount Samothrace': 1611
    }

    # The reference mountain to be excluded from the main comparison.
    olympus_height = mountains['Mount Olympus']
    print(f"For reference, the height of Mount Olympus is {olympus_height} meters.\n")
    
    # Remove Mount Olympus to find the next tallest.
    del mountains['Mount Olympus']

    # Sort the remaining mountains by elevation in descending order.
    # The 'sorted' function on dictionary items returns a list of (key, value) tuples.
    # We sort by the value (x[1]) in reverse order.
    sorted_mountains = sorted(mountains.items(), key=lambda x: x[1], reverse=True)

    # Get the name and height of the tallest mountain after Olympus.
    tallest_mountain_name = sorted_mountains[0][0]
    
    print("The historical mountains mentioned in the Iliad, ranked by height (excluding Olympus), are:")

    # Build and print the comparison string as an "equation".
    comparison_parts = []
    for name, height in sorted_mountains:
        comparison_parts.append(f"{name} ({height} m)")
    
    equation_string = " > ".join(comparison_parts)
    print(equation_string)

    print(f"\nTherefore, after Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")

find_tallest_mountain()