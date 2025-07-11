import math

def calculate_average_height():
    """
    Calculates the average adult height of the population based on genetic and nutritional information.
    """
    # Step 1: Define heights and convert them to inches.
    print("Step 1: Convert heights from feet and inches to inches.")
    height_with_milk_ft, height_with_milk_in = 4, 6
    height_without_milk_ft, height_without_milk_in = 3, 6

    height_with_milk_total_inches = height_with_milk_ft * 12 + height_with_milk_in
    height_without_milk_total_inches = height_without_milk_ft * 12 + height_without_milk_in

    print(f"Height of an individual with milk: 4'6\" = {height_with_milk_total_inches} inches.")
    print(f"Height of an individual without milk: 3'6\" = {height_without_milk_total_inches} inches.\n")

    # Step 2: Determine the proportion of the population with and without milk.
    print("Step 2: Determine the proportion of the population that grows up with versus without milk.")
    # Milk availability depends on the father's genotype. If the father is 0/0, the family has no milk.
    # The frequency of the 0/0 genotype is given as 0.5.
    freq_0_0_genotype = 0.5
    
    # The proportion of the population that grows up without milk is equal to the frequency of 0/0 fathers.
    prop_without_milk = freq_0_0_genotype
    # The rest of the population grows up with milk.
    prop_with_milk = 1 - prop_without_milk

    print(f"The proportion of fathers with genotype 0/0 (who cannot provide milk) is {prop_without_milk}.")
    print(f"Therefore, the proportion of the population growing up without milk is {prop_without_milk}.")
    print(f"The proportion of the population growing up with milk is 1 - {prop_without_milk} = {prop_with_milk}.\n")

    # Step 3: Calculate the weighted average height of the population.
    print("Step 3: Calculate the weighted average height.")
    print("The average height is the sum of each group's height multiplied by its proportion in the population.")
    print("Equation: (Proportion with milk * Height with milk) + (Proportion without milk * Height without milk)")
    
    # Calculate the average height
    avg_height = (prop_with_milk * height_with_milk_total_inches) + (prop_without_milk * height_without_milk_total_inches)
    
    # Print the equation with the numbers plugged in
    term1 = prop_with_milk * height_with_milk_total_inches
    term2 = prop_without_milk * height_without_milk_total_inches
    print(f"Average Height = ({prop_with_milk} * {height_with_milk_total_inches}) + ({prop_without_milk} * {height_without_milk_total_inches})")
    print(f"Average Height = {term1} + {term2}")
    print(f"Average Height = {avg_height} inches.\n")
    
    # Step 4: Format the result to four significant figures.
    # The result is 48. To express this with four significant figures, we write 48.00.
    final_answer_str = f"{avg_height:.2f}"
    print(f"The final answer rounded to four significant figures is {final_answer_str} inches.")
    
    # Do not remove the '<<<' and '>>>' markers.
    # They are used for extracting the final answer.
    print(f"\n<<<{final_answer_str}>>>")

# Run the calculation
calculate_average_height()