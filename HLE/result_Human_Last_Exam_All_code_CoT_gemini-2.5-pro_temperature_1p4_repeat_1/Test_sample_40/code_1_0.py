def calculate_average_height():
    """
    Calculates the average adult height of the population based on genetic and nutritional information.
    """
    # Step 1: Define heights in inches
    # Height with milk: 4 ft 6 in = 54 inches
    height_with_milk = 4 * 12 + 6
    
    # Height without milk: 3 ft 6 in = 42 inches
    height_without_milk = 3 * 12 + 6
    
    # Step 2: Determine proportions
    # Frequency of the 0/0 genotype (cannot provide milk) is given as 0.5
    proportion_without_milk = 0.5
    
    # The proportion with milk is 1 minus the proportion without milk
    proportion_with_milk = 1.0 - proportion_without_milk
    
    # Step 3: Calculate the average height
    average_height = (proportion_with_milk * height_with_milk) + (proportion_without_milk * height_without_milk)
    
    # The result is 48. To present it with four significant figures, we format it to two decimal places (48.00).
    formatted_average_height = f"{average_height:.2f}"
    
    # Print the explanation and the calculation
    print("The average adult height is calculated as a weighted average.")
    print(f"Height with milk: {height_with_milk} inches.")
    print(f"Height without milk: {height_without_milk} inches.")
    print(f"Proportion of the population that receives milk: {proportion_with_milk}")
    print(f"Proportion of the population that does not receive milk: {proportion_without_milk}")
    print("\nFinal Equation:")
    print(f"({proportion_with_milk} * {height_with_milk}) + ({proportion_without_milk} * {height_without_milk}) = {formatted_average_height}")
    
    print(f"\nThe average adult height of the population is {formatted_average_height} inches.")

# Run the calculation
calculate_average_height()
