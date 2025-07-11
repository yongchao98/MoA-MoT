def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    found in the common allotropes of pure phosphorus.
    """
    # Dictionary mapping color categories to a count of 1 for summation
    colors = {
        "White/Yellow": 1,
        "Red": 1,
        "Violet": 1,
        "Black": 1
    }

    # Get the list of numbers for the equation
    equation_numbers = list(colors.values())
    
    # Calculate the total number of colors
    total_colors = sum(equation_numbers)
    
    # Create the string for the equation
    equation_str = " + ".join(map(str, equation_numbers))

    print(f"The primary colors observed in pure phosphorus allotropes are White/Yellow, Red, Violet, and Black.")
    print(f"The total number of distinct colors is calculated as: {equation_str} = {total_colors}")

count_phosphorus_colors()