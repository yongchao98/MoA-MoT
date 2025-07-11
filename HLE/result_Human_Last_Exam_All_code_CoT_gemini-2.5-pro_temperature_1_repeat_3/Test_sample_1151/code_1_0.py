import math

def solve_pearl_riddle():
    """
    This function solves the pearl necklace riddle step-by-step.
    """
    # Part 1: Find the total number of pearls

    # First, calculate the number of pearls remaining on the string.
    # "a seven shy of eleven times eleven"
    pearls_on_string = (11 * 11) - 7

    # The fractions of the total pearls (x) that fell are:
    # 1/6, 1/5, 1/3, and 1/10.
    # To solve for the total pearls 'x', we set up the equation:
    # x = pearls_on_string + (1/6)x + (1/5)x + (1/3)x + (1/10)x
    # Rearranging the equation to solve for x:
    # x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = pearls_on_string
    # x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = pearls_on_string
    # The common denominator for 6, 5, 3, 10 is 30.
    # 1/6 -> 5/30
    # 1/5 -> 6/30
    # 1/3 -> 10/30
    # 1/10 -> 3/30
    # Sum of numerators = 5 + 6 + 10 + 3 = 24
    # So, the sum of fractions is 24/30.
    # x * (1 - 24/30) = pearls_on_string
    # x * (6/30) = pearls_on_string
    # x * (1/5) = pearls_on_string
    # x = pearls_on_string * 5
    total_pearls = pearls_on_string * 5

    print("--- Part 1: How many pearls were there altogether? ---\n")
    print("First, let's determine the number of pearls remaining on the string from the line:")
    print("'a seven shy of eleven times eleven'.")
    print(f"Calculation: (11 * 11) - 7 = 121 - 7 = {pearls_on_string}\n")

    print("Let 'x' be the total number of pearls. The full equation is the sum of the pearls on the string and the fallen parts:")
    print(f"x = {pearls_on_string} + (1/6)x + (1/5)x + (1/3)x + (1/10)x\n")
    print("To solve for x, we rearrange the equation:")
    print(f"x - (1/6 + 1/5 + 1/3 + 1/10)x = {pearls_on_string}")
    print("Finding a common denominator (30) for the fractions gives:")
    print(f"x * (1 - (5/30 + 6/30 + 10/30 + 3/30)) = {pearls_on_string}")
    print(f"x * (1 - 24/30) = {pearls_on_string}")
    print(f"x * (6/30) = {pearls_on_string}")
    print(f"x * (1/5) = {pearls_on_string}")
    print(f"x = {pearls_on_string} * 5")
    print(f"x = {total_pearls}\n")
    print(f"Therefore, there were {total_pearls} pearls altogether.\n")

    print("--------------------------------------------------\n")

    # Part 2: Find how many more pearls are needed
    fallen_pearls = total_pearls - pearls_on_string
    found_pearls = math.floor(fallen_pearls / 3)
    needed_pearls = fallen_pearls - found_pearls

    print("--- Part 2: How many more pearls are they gonna need? ---\n")
    print("First, calculate the total number of pearls that fell off the string:")
    print(f"Fallen pearls = Total pearls - Pearls on string")
    print(f"Fallen pearls = {total_pearls} - {pearls_on_string} = {fallen_pearls}\n")

    print("They find back 1/3rd of the fallen pearls:")
    print(f"Found pearls = 1/3 * {fallen_pearls} = {found_pearls}\n")

    print("The number of pearls they still need is the number of fallen pearls minus the found pearls:")
    print(f"Needed pearls = {fallen_pearls} - {found_pearls} = {needed_pearls}\n")
    print(f"Therefore, they will need {needed_pearls} more pearls.")

solve_pearl_riddle()
<<<Total pearls: 570, Pearls needed: 304>>>