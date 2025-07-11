import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle based on the provided story.
    """
    
    # Part 1: How many pearls were there altogether?
    
    # Step 1: Calculate the number of pearls that remained on the string.
    # "seven shy of eleven times eleven"
    pearls_on_string = (11 * 11) - 7
    
    # Step 2: Formulate and solve the equation for the total number of pearls (x).
    # The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string
    # To solve for x, we rearrange it:
    # x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = pearls_on_string
    # x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = pearls_on_string
    # The sum of fractions 1/6 + 1/5 + 1/3 + 1/10 = 5/30 + 6/30 + 10/30 + 3/30 = 24/30 = 4/5
    fraction_fallen = 1/6 + 1/5 + 1/3 + 1/10
    
    # x * (1 - 4/5) = pearls_on_string
    # x * (1/5) = pearls_on_string
    # x = pearls_on_string * 5
    total_pearls = pearls_on_string / (1 - fraction_fallen)
    
    # Ensure the result is an integer
    total_pearls = int(total_pearls)

    # Calculate the number of pearls for each part mentioned in the riddle
    floor_pearls = int(total_pearls * (1/6))
    bed_pearls = int(total_pearls * (1/5))
    woman_pearls = int(total_pearls * (1/3))
    lover_pearls = int(total_pearls * (1/10))

    print("--- Part 1: The Original Necklace ---")
    print(f"First, let's solve the equation of all the parts of the necklace:")
    print(f"{floor_pearls} (on the floor) + {bed_pearls} (on the bed) + {woman_pearls} (saved by woman) + {lover_pearls} (caught by lover) + {pearls_on_string} (on the string) = {total_pearls}")
    print(f"\nThere were {total_pearls} pearls altogether.\n")

    # Part 2: How many more are needed?
    
    # Step 1: Calculate the total number of fallen pearls
    fallen_pearls = total_pearls - pearls_on_string
    
    # Step 2: Calculate how many they found (1/3rd of the fallen ones)
    found_pearls = math.floor(fallen_pearls / 3)
    
    # Step 3: Calculate how many are still missing to complete the necklace
    pearls_needed = total_pearls - (pearls_on_string + found_pearls)
    
    print("--- Part 2: Rebuilding the Necklace ---")
    print(f"Total fallen pearls: {fallen_pearls}")
    print(f"Pearls they found back (1/3 of fallen): {found_pearls}")
    print(f"They now have {pearls_on_string + found_pearls} pearls in total ({pearls_on_string} on the string + {found_pearls} found).")
    print(f"\nTo restore the necklace to its original {total_pearls} pearls, they will need {pearls_needed} more.")

solve_pearl_riddle()