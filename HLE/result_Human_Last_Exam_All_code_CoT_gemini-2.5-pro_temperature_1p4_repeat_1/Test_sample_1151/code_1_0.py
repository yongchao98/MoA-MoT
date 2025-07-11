import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle step-by-step.
    """
    # --- Part 1: How many pearls altogether? ---

    print("Step 1: Determine the number of pearls remaining on the string.")
    remaining_on_string = 11 * 11 - 7
    print(f"The number of pearls remaining on the string is (11 * 11) - 7 = {remaining_on_string}.\n")

    print("Step 2: Determine the total number of pearls on the necklace.")
    print("Let X be the total number of pearls. The riddle states the fallen pearls are fractions of X:")
    print(" - 1/6 fell to the floor")
    print(" - 1/5 fell on the bed")
    print(" - 1/3 were saved by the woman")
    print(" - 1/10 were caught by her lover")
    
    print("\nThe equation for the total pearls (X) is:")
    print(f"X = (1/6 * X) + (1/5 * X) + (1/3 * X) + (1/10 * X) + {remaining_on_string}")
    
    # The sum of fractions is 1/6 + 1/5 + 1/3 + 1/10 = 5/30 + 6/30 + 10/30 + 3/30 = 24/30 = 4/5
    print("Combining the fractions gives: X = (4/5 * X) + 114")
    print("Solving for X:")
    print("X - (4/5 * X) = 114")
    print("1/5 * X = 114")
    print("X = 114 * 5")
    
    # Calculate total pearls based on the derived equation
    total_pearls = remaining_on_string * 5
    print(f"So, the total number of pearls altogether was {total_pearls}.\n")

    print("Step 3: Verify the total by showing the final equation with each part's value.")
    p_floor = int(total_pearls / 6)
    p_bed = int(total_pearls / 5)
    p_woman = int(total_pearls / 3)
    p_lover = int(total_pearls / 10)
    
    print("Let's check the numbers:")
    print(f"  {p_floor:3} pearls on the floor (1/6)")
    print(f"  {p_bed:3} pearls on the bed (1/5)")
    print(f"  {p_woman:3} pearls with the woman (1/3)")
    print(f"  {p_lover:3} pearls with the lover (1/10)")
    print(f"+ {remaining_on_string:3} pearls remaining on the string")
    print("  ---------------------------------")
    # This is the final equation with each number outputted
    print(f"  {total_pearls} = {p_floor} + {p_bed} + {p_woman} + {p_lover} + {remaining_on_string}")
    print(f"The sum of the parts is {p_floor + p_bed + p_woman + p_lover + remaining_on_string}, which matches the total.\n")

    # --- Part 2: How many more pearls do they need? ---

    print("Step 4: Calculate how many more pearls they need.")
    fallen_pearls = total_pearls - remaining_on_string
    print(f"The number of fallen pearls is {total_pearls} - {remaining_on_string} = {fallen_pearls}.")
    
    found_pearls = math.floor(fallen_pearls / 3)
    print(f"They find 1/3 of the fallen pearls: {fallen_pearls} / 3 = {found_pearls}.")
    
    needed_pearls = fallen_pearls - found_pearls
    print(f"The number of pearls they still need to find or replace is the number of fallen pearls ({fallen_pearls}) minus the number they found ({found_pearls}).")
    print(f"Therefore, they need {needed_pearls} more pearls.")
    
    # Final answer in the specified format
    print(f"\n<<< {needed_pearls} >>>")

solve_pearl_riddle()