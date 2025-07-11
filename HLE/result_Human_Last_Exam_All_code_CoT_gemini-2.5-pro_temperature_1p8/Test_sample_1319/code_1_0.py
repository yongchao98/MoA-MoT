def find_tallest_mountain_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, and prints a comparison.
    """
    # Step 1 & 2: A dictionary of historical mountains from the Iliad and their elevations in meters.
    mountains_in_iliad = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida (Turkey)": 1774,
        "Mount Samothrace (Fengari)": 1611
    }

    print("--- Mountains mentioned in the Iliad by height ---")
    print(f"The highest peak is Mount Olympus at {mountains_in_iliad['Mount Olympus']} meters.")
    
    # Step 3: Remove Mount Olympus to find the next tallest.
    mountains_without_olympus = {name: height for name, height in mountains_in_iliad.items() if name != "Mount Olympus"}

    # Step 4: Sort the remaining mountains by height in descending order.
    # The key=lambda item: item[1] tells sorted() to use the second value of each item (the height) for sorting.
    sorted_mountains = sorted(mountains_without_olympus.items(), key=lambda item: item[1], reverse=True)

    print("\nComparing the other historical mountains:")
    
    # Step 5: Construct and print the final equation/comparison.
    # This loop prints each mountain and its height.
    # The `>` character is added between them to form a visual equation.
    equation_parts = []
    for name, height in sorted_mountains:
        equation_parts.append(f"{name} ({height}m)")
    
    final_equation = " > ".join(equation_parts)
    print("Final Equation: " + final_equation)
    
    # Identify and state the tallest mountain from the remaining list.
    if sorted_mountains:
        tallest_mountain = sorted_mountains[0]
        print(f"\nAfter Mount Olympus, the tallest historical mountain is {tallest_mountain[0]}.")
    else:
        print("\nNo other mountains to compare.")

find_tallest_mountain_iliad()