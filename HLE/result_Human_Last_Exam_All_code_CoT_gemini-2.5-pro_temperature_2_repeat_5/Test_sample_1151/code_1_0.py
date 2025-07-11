def solve_pearl_riddle():
    """
    This function solves the pearl necklace riddle based on the provided text.
    """
    # Step 1 & 2: Determine the number of pearls remaining on the string.
    # "a seven shy of eleven times eleven pearls remained upon the string"
    # 11 * 11 = 121
    # 121 - 7 = 114
    remained_on_string = (11 * 11) - 7

    # Step 3: Formulate and solve the equation for the total number of pearls (x).
    # The equation is: x = (x/6) + (x/5) + (x/3) + (x/10) + 114
    # To solve, we move all 'x' terms to one side:
    # x - x/6 - x/5 - x/3 - x/10 = 114
    # To combine the fractions, we find a common denominator, which is 30.
    # (30x/30) - (5x/30) - (6x/30) - (10x/30) - (3x/30) = 114
    # (30 - 5 - 6 - 10 - 3)x / 30 = 114
    # 6x / 30 = 114
    # x / 5 = 114
    # x = 114 * 5
    total_pearls = remained_on_string * 5

    print(f"--- Solving the Riddle ---")
    print(f"\nFirst, we calculate the pearls left on the string:")
    print(f"(11 * 11) - 7 = {remained_on_string} pearls.")

    # Step 4: Display the solved equation and the total number of pearls.
    # Calculate each fractional part based on the total.
    floor_pearls = total_pearls / 6
    bed_pearls = total_pearls / 5
    woman_saved = total_pearls / 3
    lover_caught = total_pearls / 10

    print(f"\nThe equation to find the total pearls is based on adding all the parts:")
    print(f"Total = (Fell to floor) + (Fell on bed) + (Woman saved) + (Lover caught) + (Remained on string)")
    print(f"The final equation with the numbers is:")
    print(f"{total_pearls} = {int(floor_pearls)} + {int(bed_pearls)} + {int(woman_saved)} + {int(lover_caught)} + {int(remained_on_string)}")
    
    print(f"\nQuestion 1: How many pearls were there altogether?")
    print(f"Answer: There were {total_pearls} pearls altogether on the necklace.")
    
    # Step 5: Calculate how many more pearls they need.
    fallen_pearls = total_pearls - remained_on_string
    found_pearls = fallen_pearls / 3
    needed_pearls = fallen_pearls - found_pearls
    
    print(f"\nQuestion 2: How many more are they now gonna need if they manage to find back only 1/3rd of the fallen ones?")
    print(f"Total fallen pearls: {int(fallen_pearls)}")
    print(f"Pearls they found (1/3 of fallen): {int(found_pearls)}")
    print(f"Answer: They will need {int(needed_pearls)} more pearls.")

solve_pearl_riddle()
<<<304>>>