def calculate_average_height():
    """
    Calculates the average adult height of a fictional species based on genetic and environmental factors.
    """
    # Step 1: Define the proportions of the population with and without milk.
    # The father's genotype determines milk availability.
    # Frequency of fathers with 0/0 genotype (cannot provide milk) is given as 0.5.
    freq_without_milk = 0.5
    
    # The remaining proportion of fathers can provide milk.
    freq_with_milk = 1 - freq_without_milk

    # Step 2: Define the heights in inches.
    # Height with milk: 4 feet 6 inches
    height_with_milk_inches = 4 * 12 + 6
    
    # Height without milk: 3 feet 6 inches
    height_without_milk_inches = 3 * 12 + 6

    # Step 3: Calculate the average height as a weighted average.
    avg_height = (freq_with_milk * height_with_milk_inches) + (freq_without_milk * height_without_milk_inches)

    # Step 4: Print the equation and the final result.
    # The problem requires printing the final equation with the numbers.
    print("The average height is calculated as a weighted average of the heights of individuals with and without milk during their upbringing.")
    print("The proportion of the population with milk is determined by the frequency of fathers who can digest lactose.")
    print(f"Proportion with milk = {freq_with_milk}")
    print(f"Proportion without milk = {freq_without_milk}")
    print(f"Height with milk = {height_with_milk_inches} inches")
    print(f"Height without milk = {height_without_milk_inches} inches")
    print("\nFinal Equation:")
    print(f"Average Height = ({freq_with_milk} * {height_with_milk_inches}) + ({freq_without_milk} * {height_without_milk_inches})")
    
    # The result 48 needs to be formatted to four significant figures, which is 48.00.
    print(f"\nCalculated Average Height: {avg_height:.2f} inches")

calculate_average_height()
<<<48.00>>>