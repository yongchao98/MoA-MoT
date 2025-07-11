import sys
from fractions import Fraction

# Redirect print to a string to format it nicely for the final output
# This is not strictly necessary for the calculation but helps control the output format.

def solve_pearl_riddle():
    """
    This function solves the pearl riddle based on the provided text.
    It first calculates the total number of pearls on the necklace
    and then calculates how many more are needed based on the second question.
    """
    
    # --- Part 1: Calculate the total number of pearls ---
    
    # Calculate the number of pearls that remained on the string
    # "a seven shy of eleven times eleven"
    remaining_on_string = 11 * 11 - 7

    # The fractions of the total pearls that fell
    fractions_fallen = [Fraction(1, 6), Fraction(1, 5), Fraction(1, 3), Fraction(1, 10)]
    
    # The equation for the total number of pearls 'x' is:
    # x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + (11*11 - 7)
    
    print("--- Part 1: How many pearls were there altogether? ---")
    print("\nFirst, we set up the equation based on the riddle. Let 'x' be the total number of pearls.")
    print(f"The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + ({11}*11 - {7})")

    # Calculate the total fraction of pearls that fell
    total_fraction_fallen = sum(fractions_fallen)
    
    # The fraction remaining is 1 - total_fraction_fallen
    fraction_remaining = 1 - total_fraction_fallen
    
    # The number of pearls remaining on the string is equal to the total number of pearls
    # multiplied by the fraction of pearls that remained.
    # So, remaining_on_string = x * fraction_remaining
    # Which means x = remaining_on_string / fraction_remaining
    total_pearls = remaining_on_string / fraction_remaining

    print(f"\nStep 1: Calculate the number of pearls remaining on the string.")
    print(f"11 * 11 - 7 = 121 - 7 = {remaining_on_string}")

    print(f"\nStep 2: Sum the fractions of fallen pearls.")
    print(f"1/6 + 1/5 + 1/3 + 1/10 = {total_fraction_fallen}")
    
    print("\nStep 3: Rearrange the equation to solve for x.")
    print(f"x = ({total_fraction_fallen})x + {remaining_on_string}")
    print(f"x - ({total_fraction_fallen})x = {remaining_on_string}")
    print(f"({fraction_remaining})x = {remaining_on_string}")
    print(f"x = {remaining_on_string} / {fraction_remaining}")
    
    total_pearls = int(total_pearls)
    print(f"\nx = {total_pearls}")
    print(f"\nAnswer 1: There were {total_pearls} pearls on the necklace altogether.")
    
    # --- Part 2: Calculate how many more pearls are needed ---
    
    print("\n\n--- Part 2: How many more pearls are needed? ---")
    
    # Number of pearls that fell from the string
    fallen_pearls = total_pearls - remaining_on_string
    print(f"\nStep 1: Calculate the total number of fallen pearls.")
    print(f"Total Pearls - Remaining on String = {total_pearls} - {remaining_on_string} = {fallen_pearls}")

    # Number of pearls they found (1/3 of the fallen ones)
    found_pearls = fallen_pearls * Fraction(1, 3)
    found_pearls = int(found_pearls)
    print(f"\nStep 2: Calculate the number of pearls they found back (1/3 of fallen).")
    print(f"{fallen_pearls} / 3 = {found_pearls}")
    
    # Total pearls they have now
    current_total = remaining_on_string + found_pearls
    print(f"\nStep 3: Calculate the total number of pearls they currently possess.")
    print(f"Remaining on String + Found Pearls = {remaining_on_string} + {found_pearls} = {current_total}")
    
    # Number of pearls they still need
    needed_pearls = total_pearls - current_total
    print(f"\nStep 4: Calculate how many more pearls they need.")
    print(f"Original Total - Current Total = {total_pearls} - {current_total} = {needed_pearls}")
    
    print(f"\nAnswer 2: They will need {needed_pearls} more pearls.")
    
    return needed_pearls

# Execute the function and capture the final answer for the '<<<' format
final_answer = solve_pearl_riddle()
# The '<<<' part will be appended by the system based on the final returned value.
# To be compliant, let's print it directly as well.
# sys.stdout.write(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_pearl_riddle()