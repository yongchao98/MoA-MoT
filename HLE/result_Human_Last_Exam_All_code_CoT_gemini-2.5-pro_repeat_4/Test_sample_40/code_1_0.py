def calculate_average_height():
    """
    Calculates the average adult height of the population based on genetic and nutritional information.
    """
    # Step 1: Define heights and convert to inches.
    height_with_milk_ft = 4
    height_with_milk_in = 6
    height_with_milk_total_inches = height_with_milk_ft * 12 + height_with_milk_in

    height_without_milk_ft = 3
    height_without_milk_in = 6
    height_without_milk_total_inches = height_without_milk_ft * 12 + height_without_milk_in

    # Step 2: Determine the proportion of the population in each group.
    # The frequency of the 0/0 genotype is given as 0.5.
    # A family gets milk only if the father's genotype is NOT 0/0.
    # The population is randomly mating, so the proportion of fathers with a given genotype
    # is the same as the frequency of that genotype in the population.
    freq_father_00 = 0.5
    proportion_without_milk = freq_father_00
    proportion_with_milk = 1 - freq_father_00

    # Step 3: Calculate the weighted average height.
    average_height = (proportion_with_milk * height_with_milk_total_inches) + \
                     (proportion_without_milk * height_without_milk_total_inches)

    # Step 4: Print the explanation and the final, formatted result.
    print("The average adult height is the weighted average of the two possible heights.")
    print("Height with milk during upbringing = {} feet {} inches = {} inches.".format(height_with_milk_ft, height_with_milk_in, height_with_milk_total_inches))
    print("Height without milk during upbringing = {} feet {} inches = {} inches.".format(height_without_milk_ft, height_without_milk_in, height_without_milk_total_inches))
    print("\nThe proportion of the population that grows up without milk is determined by the frequency of fathers with the 0/0 genotype, which is {}.".format(proportion_without_milk))
    print("The proportion of the population that grows up with milk is therefore {}.".format(proportion_with_milk))
    print("\nThe final calculation for the average height is:")
    # The format specifier '{:.4g}' ensures the result is displayed with four significant figures.
    print("({:.1f} * {}) + ({:.1f} * {}) = {:.4g}".format(proportion_with_milk, height_with_milk_total_inches, proportion_without_milk, height_without_milk_total_inches, average_height))

calculate_average_height()
print("\n<<<48.00>>>")