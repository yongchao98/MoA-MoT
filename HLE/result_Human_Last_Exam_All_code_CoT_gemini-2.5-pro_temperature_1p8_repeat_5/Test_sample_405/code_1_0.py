def solve_phosphorus_colors():
    """
    This function identifies the colors of pure phosphorus allotropes and calculates the total number.
    """
    # Step 1 & 2: A dictionary storing the major allotropes of phosphorus and their colors.
    # Yellow phosphorus is often considered an impure or light-degraded form of white phosphorus,
    # so we will count white as the primary color for that form.
    allotropes_info = {
        "White Phosphorus": "white",
        "Red Phosphorus": "red",
        "Violet Phosphorus": "violet",
        "Black Phosphorus": "black"
    }

    print("The colors observed in pure allotropes of phosphorus are:")
    for allotrope, color in allotropes_info.items():
        print(f"- {color.capitalize()} (from {allotrope})")

    # Step 3: Count the number of unique colors.
    num_colors = len(allotropes_info)

    # Step 4: Create and display the final equation.
    print("\nTo find the total, we sum one for each distinct color:")

    # Create a list like ["1 (white)", "1 (red)", ...]
    equation_components = [f"1 ({color})" for color in allotropes_info.values()]

    # Join the components with " + " to form the equation string
    equation_string = " + ".join(equation_components)

    # Print the full equation showing each number
    print(f"{equation_string} = {num_colors}")

solve_phosphorus_colors()