import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle based on the provided text.
    """
    # Part 1: How many pearls were there altogether?

    # First, calculate the pearls remaining on the string:
    # "seven shy of eleven times eleven"
    remaining_on_string = 11 * 11 - 7

    # Let the total number of pearls be 'x'. The riddle states:
    # x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_on_string
    # To solve for x, we first find the sum of the fractions:
    # 1/6 + 1/5 + 1/3 + 1/10 = 5/30 + 6/30 + 10/30 + 3/30 = 24/30 = 4/5
    # So, the equation is: x = (4/5)x + remaining_on_string
    # x - (4/5)x = remaining_on_string
    # (1/5)x = remaining_on_string
    # x = remaining_on_string * 5
    
    total_pearls = remaining_on_string * 5

    # Calculate the number of pearls for each part for the final equation printout
    floor_pearls = total_pearls / 6
    bed_pearls = total_pearls / 5
    woman_pearls = total_pearls / 3
    lover_pearls = total_pearls / 10

    print("--- Part 1: Total Pearls ---")
    print("The equation to find the total number of pearls (x) is based on the fractions given:")
    # Using the calculated numbers to show the final solved equation
    print(f"Total Pearls ({int(total_pearls)}) = On the Floor ({int(floor_pearls)}) + On the Bed ({int(bed_pearls)}) + Saved by Woman ({int(woman_pearls)}) + Caught by Lover ({int(lover_pearls)}) + Remained on String ({int(remaining_on_string)})")
    
    print(f"\nThere were {int(total_pearls)} pearls altogether.\n")


    # Part 2: How many more are they now gonna need?

    # Calculate the total number of pearls that fell from the string
    total_fallen = total_pearls - remaining_on_string

    # Calculate how many of the fallen pearls they found (1/3)
    found_pearls = total_fallen / 3

    # The number they still need is the number of pearls that are still missing
    needed_pearls = total_fallen - found_pearls
    
    print("--- Part 2: Needed Pearls ---")
    print(f"Total number of fallen pearls: {int(total_fallen)}")
    print(f"Number of pearls they found back (1/3 of fallen): {int(found_pearls)}")
    print(f"Therefore, the number of pearls they still need to find is {int(needed_pearls)}.")


# Execute the function to solve the riddle
solve_pearl_riddle()