def find_tallest_mountain():
    """
    Finds the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Tmolus": 2157,
        "Mount Athos": 2033,
        "Mount Ossa": 1978,
        "Mount Ida (Trojan)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    # Remove Mount Olympus from the comparison as per the instructions.
    del mountains["Mount Olympus"]

    # Find the name of the tallest remaining mountain.
    tallest_mountain_name = max(mountains, key=mountains.get)
    
    # Sort the mountains by height in descending order for a clear comparison.
    sorted_mountains = sorted(mountains.items(), key=lambda item: item[1], reverse=True)
    
    print("The tallest historical mountain in the Iliad (after Olympus) is found by comparing the elevations:")
    
    # Constructing the comparison equation string.
    equation_parts = []
    for name, height in sorted_mountains:
        print(f"- {name}: {height}m")
        equation_parts.append(str(height))
    
    # The final equation showing the comparison of all heights.
    comparison_equation = " > ".join(equation_parts)
    print(f"\nFinal Comparison Equation:")
    print(comparison_equation)

    print(f"\nTherefore, the tallest mountain is {tallest_mountain_name}.")

find_tallest_mountain()