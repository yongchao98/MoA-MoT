def calculate_average_height():
    """
    Calculates the average adult height of a fictional mammalian species
    based on genetic and environmental factors.
    """
    # 1. Define variables based on the problem description.
    # The frequency of the 0/0 genotype determines the proportion of fathers
    # who cannot provide milk. This is given as half the population.
    freq_father_is_00 = 0.5

    # This directly translates to the proportion of the population raised with or without milk.
    prop_without_milk = freq_father_is_00
    prop_with_milk = 1 - prop_without_milk

    # Heights with and without milk are given.
    # With milk: 1 foot taller than without milk.
    # Without milk: 3 feet and 6 inches.
    height_without_milk_ft = 3
    height_without_milk_in_part = 6
    # With milk: 3 feet 6 inches + 1 foot = 4 feet 6 inches
    height_with_milk_ft = 4
    height_with_milk_in_part = 6

    # 2. Convert all height measurements to inches for calculation.
    # (1 foot = 12 inches)
    height_with_milk_inches = height_with_milk_ft * 12 + height_with_milk_in_part
    height_without_milk_inches = height_without_milk_ft * 12 + height_without_milk_in_part

    # 3. Calculate the weighted average height for the entire population.
    part1_calculation = prop_with_milk * height_with_milk_inches
    part2_calculation = prop_without_milk * height_without_milk_inches
    average_height = part1_calculation + part2_calculation

    # 4. Print the step-by-step logic and the final calculation.
    print("The average adult height depends on whether an individual had milk during their upbringing.")
    print("This is determined by their father's genotype, as he is the sole milk provider.")
    print(f"\nThe frequency of fathers with the 0/0 genotype (unable to provide milk) is {freq_father_is_00}.")
    print(f"Therefore, the proportion of the population raised WITHOUT milk is {prop_without_milk}.")
    print(f"The proportion of the population raised WITH milk is {prop_with_milk}.")
    print("\n---")
    print("Height breakdown:")
    print(f"Height with milk: {height_with_milk_ft} ft {height_with_milk_in_part} in = {height_with_milk_inches} inches.")
    print(f"Height without milk: {height_without_milk_ft} ft {height_without_milk_in_part} in = {height_without_milk_inches} inches.")
    print("---\n")
    print("The average population height is the weighted average:")
    print(f"Average Height = (Proportion with milk * Height with milk) + (Proportion without milk * Height without milk)")
    print(f"Average Height = ({prop_with_milk} * {height_with_milk_inches}) + ({prop_without_milk} * {height_without_milk_inches})")
    print(f"Average Height = {part1_calculation} + {part2_calculation}")
    
    # Final answer formatted to four significant figures.
    print(f"Final Average Height = {average_height:.2f} inches.")

calculate_average_height()
<<<48.00>>>