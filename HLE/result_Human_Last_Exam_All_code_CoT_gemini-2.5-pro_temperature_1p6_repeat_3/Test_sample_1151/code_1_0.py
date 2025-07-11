import math

def solve_pearl_riddle():
    """
    This function solves the pearl necklace riddle by setting up and solving
    a linear equation. It then calculates the number of pearls needed for replacement.
    """
    # Part 1: How many pearls were there altogether?

    # First, calculate the number of pearls that remained on the string.
    # The riddle states "a seven shy of eleven times eleven".
    eleven = 11
    seven = 7
    remaining_pearls = (eleven * eleven) - seven

    # The equation for the total number of pearls (x) is:
    # x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_pearls
    # To solve for x, we rearrange the equation:
    # x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = remaining_pearls
    # x * (1 - 1/6 - 1/5 - 1/3 - 1/10) = remaining_pearls
    # The common denominator for 6, 5, 3, and 10 is 30.
    # x * (30/30 - 5/30 - 6/30 - 10/30 - 3/30) = remaining_pearls
    # x * (6/30) = remaining_pearls
    # x * (1/5) = remaining_pearls
    # x = remaining_pearls * 5

    total_pearls = remaining_pearls * 5

    # Now we calculate the number of pearls in each fallen group for the final output.
    floor_pearls = total_pearls / 6
    bed_pearls = total_pearls / 5
    woman_saved = total_pearls / 3
    lover_caught = total_pearls / 10

    print("--- Part 1: Finding the Total Number of Pearls ---")
    print(f"The number of pearls remaining on the string is: {eleven} * {eleven} - {seven} = {remaining_pearls}")
    print("\nThe full equation for the total pearls is:")
    print(f"Total = (Fell to floor) + (Fell on bed) + (Woman saved) + (Lover caught) + (Remained on string)")
    # Using math.ceil to handle any potential floating point inaccuracies and ensure integer output
    print(f"{math.ceil(total_pearls)} = {math.ceil(floor_pearls)} + {math.ceil(bed_pearls)} + {math.ceil(woman_saved)} + {math.ceil(lover_caught)} + {math.ceil(remaining_pearls)}")
    print(f"\nThere were {math.ceil(total_pearls)} pearls altogether.")

    # Part 2: How many more pearls are needed?
    total_fallen_pearls = total_pearls - remaining_pearls
    pearls_found = total_fallen_pearls / 3
    
    # They need to replace the ones they did not find (2/3 of the fallen ones).
    pearls_needed = total_fallen_pearls - pearls_found

    print("\n--- Part 2: Calculating the Number of Pearls to Replace ---")
    print(f"Total number of fallen pearls: {math.ceil(total_fallen_pearls)}")
    print(f"They find 1/3 of the fallen pearls, so they still need to replace the other 2/3.")
    print(f"Number of pearls to buy = {math.ceil(total_fallen_pearls)} - (1/3 * {math.ceil(total_fallen_pearls)}) = {math.ceil(pearls_needed)}")
    print(f"\nThey will need {math.ceil(pearls_needed)} more pearls.")

# Run the solver
solve_pearl_riddle()

# Final answer format
print("\n<<<570, 304>>>")