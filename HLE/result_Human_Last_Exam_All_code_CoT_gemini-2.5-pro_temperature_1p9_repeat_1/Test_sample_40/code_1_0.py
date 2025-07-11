def calculate_average_height():
    """
    Calculates the average adult height of a fictional mammalian population
    based on genetic and nutritional factors.
    """
    # Step 1: Define the given proportions and heights.
    # The proportion of the population with genotype 0/0 is 0.5.
    # A child receives milk only if their father's genotype is not 0/0.
    # Therefore, the proportion of the population receiving milk is 1 - 0.5.
    prop_without_milk = 0.5
    prop_with_milk = 1 - prop_without_milk

    # Heights in inches.
    # Height with milk = 4 ft 6 in = (4 * 12) + 6
    height_with_milk = 54
    # Height without milk = 3 ft 6 in = (3 * 12) + 6
    height_without_milk = 42

    # Step 2: Calculate the weighted average height.
    average_height = (prop_with_milk * height_with_milk) + (prop_without_milk * height_without_milk)

    # Step 3: Print the explanation and the final result.
    print("To find the average adult height, we calculate a weighted average based on the proportion of the population that receives milk.")
    print("\nA child's nutrition depends on the father's ability to digest lactose.")
    print(f"The frequency of fathers who cannot digest milk (genotype 0/0) is given as: {prop_without_milk}")
    print(f"Therefore, the proportion of the population that grows up *without* milk is: {prop_without_milk}")
    print(f"The proportion of the population that grows up *with* milk is 1 - {prop_without_milk} = {prop_with_milk}\n")

    print(f"The average height for an individual with milk is {height_with_milk} inches.")
    print(f"The average height for an individual without milk is {height_without_milk} inches.\n")

    print("The equation for the average height is:")
    # The f-string formats the result to 4 significant figures (2 decimal places for this value).
    print(f"({prop_with_milk} * {height_with_milk}) + ({prop_without_milk} * {height_without_milk}) = {average_height:.2f} inches")

    # The final answer format required by the user prompt.
    print(f"\n<<<{average_height:.2f}>>>")

calculate_average_height()