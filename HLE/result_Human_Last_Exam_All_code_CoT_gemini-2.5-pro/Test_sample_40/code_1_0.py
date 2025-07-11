def calculate_average_height():
    """
    Calculates the average adult height of a fictional species based on genetic and environmental factors.
    """
    # The frequency of the 0/0 genotype in the population, as given in the problem.
    freq_father_no_milk = 0.5

    # The frequency of fathers who provide milk is the remainder of the population.
    freq_father_with_milk = 1 - freq_father_no_milk

    # The height of individuals who grew up without milk, converted to inches.
    # 3 feet 6 inches = (3 * 12) + 6 = 42 inches.
    height_no_milk_inches = (3 * 12) + 6

    # The height of individuals who grew up with milk, converted to inches.
    # 4 feet 6 inches = (4 * 12) + 6 = 54 inches.
    height_with_milk_inches = (4 * 12) + 6

    # Since a child's height is determined by their father's genotype, the proportion of
    # short and tall individuals in the population matches the frequency of fathers
    # who can or cannot provide milk.
    prop_short = freq_father_no_milk
    prop_tall = freq_father_with_milk

    # Calculate the weighted average height of the population.
    average_height = (prop_tall * height_with_milk_inches) + (prop_short * height_no_milk_inches)

    # Print the equation used for the calculation, showing each number.
    # The final answer is formatted to have two decimal places, giving it four significant figures.
    print(f"The average height is calculated as a weighted average of the two possible heights:")
    print(f"({prop_tall} * {height_with_milk_inches}) + ({prop_short} * {height_no_milk_inches}) = {average_height:.2f} inches")

calculate_average_height()