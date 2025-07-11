def calculate_average_height():
    """
    Calculates the average adult height of a fictional mammalian species
    based on genetic and environmental factors.
    """

    # Step 1: Define the given proportions and heights.
    # The proportion of the population with genotype 0/0 is 0.5.
    # A family has no milk if the father's genotype is 0/0.
    # Therefore, the proportion of the population raised without milk is 0.5.
    prop_without_milk = 0.5
    prop_with_milk = 1 - prop_without_milk

    # Heights in feet and inches.
    height_with_milk_ft, height_with_milk_in_part = 4, 6
    height_without_milk_ft, height_without_milk_in_part = 3, 6

    # Step 2: Convert heights to inches.
    height_with_milk_inches = height_with_milk_ft * 12 + height_with_milk_in_part
    height_without_milk_inches = height_without_milk_ft * 12 + height_without_milk_in_part

    # Step 3: Calculate the weighted average height.
    # Average Height = (Proportion with milk * Height with milk) + (Proportion without milk * Height without milk)
    avg_height = (prop_with_milk * height_with_milk_inches) + \
                 (prop_without_milk * height_without_milk_inches)

    # Step 4: Print the components of the final equation and the result.
    # The final result is formatted to two decimal places to achieve four significant figures.
    print(f"The calculation for the average height is a weighted average:")
    print(f"({prop_with_milk:.2f} of population * {height_with_milk_inches} inches) + ({prop_without_milk:.2f} of population * {height_without_milk_inches} inches)")

    result_part1 = prop_with_milk * height_with_milk_inches
    result_part2 = prop_without_milk * height_without_milk_inches
    
    print(f"This simplifies to:")
    print(f"{result_part1} + {result_part2} = {avg_height:.2f} inches")

    print(f"\nThe average adult height of the population is {avg_height:.2f} inches.")

    # Final answer in the required format
    print(f"\n<<<{avg_height:.2f}>>>")

calculate_average_height()