import math

def calculate_average_height():
    """
    Calculates the average adult height of the species based on genetic and nutritional information.
    """
    # Step 1: Define the heights in inches.
    height_with_milk_ft = 4
    height_with_milk_in = 6
    height_without_milk_ft = 3
    height_without_milk_in = 6

    # Convert heights to total inches for calculation.
    # Height with milk = 4 ft 6 in = (4 * 12) + 6 = 54 inches
    # Height without milk = 3 ft 6 in = (3 * 12) + 6 = 42 inches
    height_with_milk = height_with_milk_ft * 12 + height_with_milk_in
    height_without_milk = height_without_milk_ft * 12 + height_without_milk_in

    # Step 2: Determine the proportion of the population with and without milk.
    # The availability of milk depends on the father's genotype.
    # If the father is 0/0, the family has no milk.
    # The problem states the frequency of the 0/0 genotype is 0.5.
    proportion_of_fathers_0_0 = 0.5

    # The proportion of the population raised without milk is equal to the proportion of fathers
    # with the 0/0 genotype.
    proportion_without_milk = proportion_of_fathers_0_0
    
    # The remaining proportion of the population is raised with milk.
    proportion_with_milk = 1 - proportion_without_milk

    # Step 3: Calculate the weighted average height.
    # Average Height = (Proportion with milk * Height with milk) + (Proportion without milk * Height without milk)
    average_height = (proportion_with_milk * height_with_milk) + (proportion_without_milk * height_without_milk)
    
    # Step 4: Print the final equation and result.
    print("The average height calculation is based on the following equation:")
    print(f"({proportion_with_milk} * {height_with_milk}) + ({proportion_without_milk} * {height_without_milk}) = {average_height:.4g}")
    
    # Final result rounded to four significant figures.
    final_answer = f"{average_height:.4g}"
    # Per instructions, the final result can have trailing zeros to show 4 significant figures
    if '.' not in final_answer:
        final_answer += '.'
    while len(final_answer.replace('.', '')) < 4:
        final_answer += '0'

    print(f"\nThe average adult height of the population is {final_answer} inches.")
    return final_answer

final_answer = calculate_average_height()
# The final output needs to be wrapped according to instructions
# print(f"<<<{final_answer}>>>") # This would be part of the final wrapper if executed