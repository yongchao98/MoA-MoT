import math

def calculate_average_height():
    """
    Calculates the average adult height of the described species population.
    """
    # Step 1: Define the given parameters.
    # The proportion of the population with genotype 0/0 is 0.5.
    # This is the proportion of fathers who cannot provide milk.
    prop_without_milk = 0.5
    
    # The rest of the fathers can provide milk.
    prop_with_milk = 1 - prop_without_milk
    
    # Heights are given in feet and inches.
    height_no_milk_ft, height_no_milk_in = 3, 6
    height_with_milk_ft, height_with_milk_in = 4, 6 # 1 foot taller
    
    # Step 2: Convert heights to inches for a consistent unit.
    height_no_milk_inches = height_no_milk_ft * 12 + height_no_milk_in
    height_with_milk_inches = height_with_milk_ft * 12 + height_with_milk_in
    
    # Step 3: Calculate the weighted average height of the population.
    # This is the sum of each group's height multiplied by its proportion in the population.
    avg_height = (prop_without_milk * height_no_milk_inches) + (prop_with_milk * height_with_milk_inches)
    
    # Step 4: Print the calculation and the final result formatted to four significant figures.
    print("The average population height is calculated as a weighted average:")
    print("Average Height = (Proportion without milk * Height without milk) + (Proportion with milk * Height with milk)\n")
    
    print("The numbers for the final equation are:")
    print(f"Proportion without milk: {prop_without_milk}")
    print(f"Height without milk (inches): {height_no_milk_inches}")
    print(f"Proportion with milk: {prop_with_milk}")
    print(f"Height with milk (inches): {height_with_milk_inches}\n")

    print("Final Equation:")
    # The problem asks to output each number in the final equation.
    print(f"Average Height = ({prop_without_milk} * {height_no_milk_inches}) + ({prop_with_milk} * {height_with_milk_inches}) = {avg_height:.2f} inches")
    
    # The final answer must be formatted to four significant figures. 48.00 has four.
    final_answer_formatted = f"{avg_height:.2f}"
    
    print(f"\nThe average adult height of the population is {final_answer_formatted} inches.")
    
    # Output the final answer in the required format.
    print(f"<<<{final_answer_formatted}>>>")

calculate_average_height()