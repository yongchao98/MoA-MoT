def solve_pearl_riddle():
    """
    This function solves the pearl necklace riddle step-by-step.
    """
    # Part 1: Find the total number of pearls on the necklace.

    # First, calculate the number of pearls remaining on the string from the clue:
    # "a seven shy of eleven times eleven"
    remaining_on_string = (11 * 11) - 7

    # The equation for the total number of pearls (x) is:
    # x = (pearls on floor) + (pearls on bed) + (saved by woman) + (caught by lover) + (remaining on string)
    # x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_on_string
    # To solve for x, we can rearrange the equation:
    # x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = remaining_on_string
    # This is equivalent to:
    # x * (1 - 1/6 - 1/5 - 1/3 - 1/10) = remaining_on_string
    # The sum of the fractions is (5/30 + 6/30 + 10/30 + 3/30) = 24/30 = 4/5.
    # So, x * (1 - 4/5) = remaining_on_string
    # x * (1/5) = remaining_on_string
    # x = remaining_on_string * 5
    total_pearls = remaining_on_string * 5

    print("--- Part 1: How many pearls were there altogether? ---")
    print("Let 'x' be the total number of pearls.")
    print("First, we find the number of pearls left on the string:")
    print(f"Equation: (11 * 11) - 7 = 121 - 7 = {remaining_on_string}")
    print("\nThe main equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string")
    print(f"So, the equation with our numbers is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + {remaining_on_string}")
    print("\nSolving for x:")
    print(f"x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = {remaining_on_string}")
    print(f"Combining the fractions of x gives (24/30)x, which simplifies to (4/5)x.")
    print(f"x - (4/5)x = {remaining_on_string}")
    print(f"(1/5)x = {remaining_on_string}")
    print(f"x = {remaining_on_string} * 5")
    print(f"x = {total_pearls}")
    print(f"\nThere were {total_pearls} pearls on the necklace altogether.")

    # Part 2: Find how many more pearls are needed.
    fallen_pearls = total_pearls - remaining_on_string
    found_pearls = int(fallen_pearls / 3)
    current_pearls_count = remaining_on_string + found_pearls
    pearls_needed = total_pearls - current_pearls_count

    print("\n--- Part 2: How many more pearls are they gonna need? ---")
    print(f"First, find the total number of fallen pearls:")
    print(f"Equation: {total_pearls} (total) - {remaining_on_string} (on string) = {fallen_pearls} fallen pearls")
    print(f"\nThey find 1/3rd of the fallen pearls:")
    print(f"Equation: {fallen_pearls} / 3 = {found_pearls} found pearls")
    print(f"\nTotal pearls they have now:")
    print(f"Equation: {remaining_on_string} (on string) + {found_pearls} (found) = {current_pearls_count} pearls")
    print(f"\nNumber of pearls they still need to restore the necklace:")
    print(f"Equation: {total_pearls} (original total) - {current_pearls_count} (they have) = {pearls_needed} pearls")
    print(f"\nThey will need {pearls_needed} more pearls.")

    return pearls_needed

# Run the function and print the final answer in the required format.
final_answer = solve_pearl_riddle()
print(f"\n<<<{final_answer}>>>")
