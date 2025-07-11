import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle step by step.
    """
    # Part 1: How many pearls were there altogether?
    print("--- Part 1: How many pearls were there altogether? ---")
    
    # Step 1: Calculate the pearls remaining on the string.
    remained_on_string = (11 * 11) - 7
    print("The number of pearls remaining on the string is 'seven shy of eleven times eleven'.")
    print(f"Calculation: (11 * 11) - 7 = 121 - 7 = {remained_on_string}\n")
    
    # Step 2 & 3: Set up and solve the equation for the total number of pearls (x).
    # The equation is: x = (x/6) + (x/5) + (x/3) + (x/10) + remained_on_string
    # To solve for x, we can find the fraction of pearls that remained on the string.
    # Fraction of fallen pearls = 1/6 + 1/5 + 1/3 + 1/10 = (5+6+10+3)/30 = 24/30 = 4/5
    # Fraction remaining on string = 1 - 4/5 = 1/5
    # So, (1/5) * x = remained_on_string
    # x = remained_on_string * 5
    
    total_pearls = remained_on_string * 5

    print("Let 'x' be the total number of pearls on the necklace.")
    print("The fractions of the total that fell are 1/6, 1/5, 1/3, and 1/10.")
    print("The equation for the total number of pearls is:")
    print(f"x = (x/6) + (x/5) + (x/3) + (x/10) + {remained_on_string}")
    print("\nThis means the fraction of pearls remaining on the string is 1 - (1/6 + 1/5 + 1/3 + 1/10) = 1/5.")
    print(f"Therefore, 1/5 of the total pearls is equal to the {remained_on_string} pearls left on the string.")
    print(f"Total pearls = {remained_on_string} * 5")
    print(f"The total number of pearls on the necklace was: {total_pearls}\n")
    
    print("="*50 + "\n")

    # Part 2: How many more pearls are needed?
    print("--- Part 2: How many more pearls are needed? ---")
    
    # Step 4: Calculate the total number of fallen pearls.
    fallen_pearls = total_pearls - remained_on_string
    print(f"Total fallen pearls = {total_pearls} (total) - {remained_on_string} (on string) = {fallen_pearls} pearls.")
    
    # Step 5: Calculate how many pearls they find.
    found_pearls = math.floor(fallen_pearls / 3)
    print(f"They find 1/3 of the fallen pearls: {fallen_pearls} / 3 = {found_pearls} pearls.")

    # Step 6: Calculate how many more pearls are needed.
    needed_pearls = fallen_pearls - found_pearls
    print(f"The number of pearls still needed is the difference between what fell and what was found.")
    print(f"Needed pearls = {fallen_pearls} - {found_pearls}")
    print(f"\nThey will need {needed_pearls} more pearls to complete the necklace.")


solve_pearl_riddle()
<<<304>>>