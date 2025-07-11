def analyze_sorting_efficiency():
    """
    This function calculates the efficiency of the cell sorting experiment
    and provides an interpretation based on the results.
    """
    # Total number of wells sorted, each intended to contain one doublet
    sorted_doublets = 95

    # Number of wells from which paired heavy and light chains were successfully amplified
    successful_pairs = 4

    # Calculate the experimental efficiency
    efficiency_percentage = (successful_pairs / sorted_doublets) * 100

    print("Step 1: Define the experimental parameters.")
    print(f"Total events sorted as doublets: {sorted_doublets}")
    print(f"Successful outcomes (paired chains recovered): {successful_pairs}")
    print("-" * 50)
    
    print("Step 2: Calculate the experimental efficiency.")
    # The final code must output each number in the final equation.
    print(f"Efficiency Equation: ({successful_pairs} / {sorted_doublets}) * 100")
    print(f"Result: The final success rate is {efficiency_percentage:.2f}%.")
    print("-" * 50)

    print("Conclusion:")
    print("The extremely low efficiency of ~4.21% supports the hypothesis that most sorted events were not true biological interactions.")
    print("This is because true antigen-specific B cells are very rare, making random cell pairings the much more common event identified as a 'doublet' by the sorter.")

analyze_sorting_efficiency()