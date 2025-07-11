import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle based on the provided text.
    """
    # Part 1: Find the total number of pearls in the necklace.

    # Calculate the number of pearls remaining on the string: "seven shy of eleven times eleven"
    pearls_on_string = 11 * 11 - 7

    # The equation for the total number of pearls (x) is:
    # x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string
    # To solve for x, we rearrange the equation:
    # x - (1/6 + 1/5 + 1/3 + 1/10)x = pearls_on_string
    # The sum of the fractions is 24/30, or 4/5.
    # x - (4/5)x = pearls_on_string
    # (1/5)x = pearls_on_string
    # x = pearls_on_string * 5
    
    total_pearls = pearls_on_string * 5

    print("Step 1: Finding the total number of pearls.")
    print("The problem can be solved with the equation where 'x' is the total number of pearls:")
    # The problem asks to output each number in the final equation.
    print(f"(1/6)x + (1/5)x + (1/3)x + (1/10)x + {pearls_on_string} = x")
    print(f"By solving the equation, we find the total number of pearls altogether was: {total_pearls}")
    print("-" * 40)

    # Part 2: How many more pearls are needed.

    # Calculate the total number of fallen pearls
    fallen_pearls = total_pearls - pearls_on_string

    # They find 1/3rd of the fallen ones
    found_pearls = fallen_pearls / 3

    # The number they still need is the difference between what fell and what was found
    needed_pearls = fallen_pearls - found_pearls
    
    print("Step 2: Calculating how many more pearls they need.")
    print(f"The number of pearls that fell off the string is: {int(fallen_pearls)}")
    print(f"They manage to find 1/3 of the fallen pearls, which is: {math.floor(found_pearls)}")
    print(f"Therefore, the number of pearls they still need to find is: {int(needed_pearls)}")

solve_pearl_riddle()