def solve_graphene_classification():
    """
    This function classifies three graphene nanoribbon band structures
    and concatenates the results.
    """

    # Classification data derived from visual analysis of the plots
    # Each tuple contains: (Edge Type, Width N, Band Type)
    classifications = [
        ('A', 8, 0),  # Plot 1: Armchair, N=8, Metallic
        ('Z', 7, 0),  # Plot 2: Zigzag, N=7, Metallic
        ('A', 7, 1)   # Plot 3: Armchair, N=7, Semiconducting
    ]

    final_string = ""
    print("Analyzing and concatenating classifications for each plot:")

    for i, (edge, width, band) in enumerate(classifications):
        plot_num = i + 1
        # Format the classification string for the current plot
        class_str = f"{edge}{width}{band}"
        
        # Print the components for clarity, showing each number
        print(f"Plot {plot_num}: Edge='{edge}', Width={width}, Band={band} -> Classification='{class_str}'")
        
        # Append to the final concatenated string
        final_string += class_str

    print("\nFinal concatenated classification string:")
    print(final_string)

solve_graphene_classification()