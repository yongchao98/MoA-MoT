def find_tallest_mountain():
    # A dictionary of historical mountains mentioned in the Iliad (excluding Olympus) and their elevations in meters.
    mountains = {
        "Mount Ida (Troad)": 1774,
        "Mount Samothrace": 1611,
        "Mount Athos": 2033,
        "Mount Ossa": 1978,
        "Mount Pelion": 1624
    }

    # Find the mountain with the maximum elevation
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    # Create the equation string
    # The user wants to see each number in the final equation.
    heights = [str(h) for h in mountains.values()]
    equation = f"max({', '.join(heights)}) = {tallest_mountain_height}"
    
    print("Finding the tallest mountain mentioned in the Iliad (after Mount Olympus):")
    for name, height in mountains.items():
        print(f"- {name}: {height}m")
    
    print("\nThe comparison equation is:")
    print(equation)
    
    print(f"\nAfter Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

find_tallest_mountain()