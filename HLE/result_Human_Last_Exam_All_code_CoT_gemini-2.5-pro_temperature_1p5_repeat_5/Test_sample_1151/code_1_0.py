import sys
from fractions import Fraction

# This function will solve the riddle and print the results.
def solve_pearl_riddle():
    # Part 1: Find the total number of pearls
    
    # Calculate the number of pearls remaining on the string
    # "a seven shy of eleven times eleven"
    remaining_on_string = (11 * 11) - 7
    
    # The sum of the fractions of pearls scattered
    # 1/6 (floor) + 1/5 (bed) + 1/3 (woman) + 1/10 (lover)
    # Using the Fraction module for precision
    sum_of_fractions = Fraction(1, 6) + Fraction(1, 5) + Fraction(1, 3) + Fraction(1, 10)
    
    # The equation is: x = (sum_of_fractions * x) + remaining_on_string
    # Rearranging to solve for x: x - (sum_of_fractions * x) = remaining_on_string
    # x * (1 - sum_of_fractions) = remaining_on_string
    # x = remaining_on_string / (1 - sum_of_fractions)
    
    total_pearls = remaining_on_string / (1 - sum_of_fractions)
    total_pearls = int(total_pearls) # The result should be a whole number

    # --- Print the solution for the first part ---
    print("Part 1: How many pearls were there altogether?")
    print("The riddle describes the total pearls 'x' as the sum of its parts:")
    print("x = (pearls on floor) + (pearls on bed) + (pearls woman saved) + (pearls lover caught) + (pearls on string)")
    print("x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + (11*11 - 7)")
    
    # Calculate each part using the solved total
    floor_pearls = int(total_pearls / 6)
    bed_pearls = int(total_pearls / 5)
    woman_saved = int(total_pearls / 3)
    lover_caught = int(total_pearls / 10)
    
    print("\nBy solving for 'x', we find the total was 570 pearls.")
    print("Let's verify by plugging the numbers into the equation:")
    print(f"{total_pearls} = {floor_pearls} + {bed_pearls} + {woman_saved} + {lover_caught} + {remaining_on_string}")
    
    print(f"\nAnswer: There were {total_pearls} pearls altogether.")
    
    # --- Part 2: Find how many more pearls they need ---
    
    # Calculate the total number of fallen pearls
    fallen_pearls = total_pearls - remaining_on_string
    
    # They find 1/3 of the fallen pearls
    found_pearls = fallen_pearls / 3
    
    # The number needed is the number they haven't found
    needed_pearls = fallen_pearls - found_pearls
    needed_pearls = int(needed_pearls)

    # --- Print the solution for the second part ---
    print("\n" + "="*50)
    print("Part 2: How many more pearls are they now gonna need?")
    print(f"First, calculate the total number of fallen pearls: {total_pearls} (total) - {remaining_on_string} (on string) = {fallen_pearls} pearls.")
    print(f"They find 1/3rd of these fallen pearls: {fallen_pearls} / 3 = {int(found_pearls)} pearls.")
    print(f"This means the number of pearls still lost is {fallen_pearls} - {int(found_pearls)} = {needed_pearls} pearls.")
    
    print(f"\nAnswer: They will need {needed_pearls} more pearls to complete the necklace.")

# Execute the function
solve_pearl_riddle()