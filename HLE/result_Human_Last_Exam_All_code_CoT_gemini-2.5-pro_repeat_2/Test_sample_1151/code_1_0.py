from fractions import Fraction

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle step-by-step.
    """
    # Part 1: Calculate the total number of pearls on the necklace.

    # The number of pearls remaining on the string is "seven shy of eleven times eleven".
    eleven = 11
    seven = 7
    remaining_pearls = eleven * eleven - seven

    # The sum of the fractions of the total pearls that fell off.
    # Using the Fraction class to avoid floating-point inaccuracies.
    fraction_fallen = Fraction(1, 6) + Fraction(1, 5) + Fraction(1, 3) + Fraction(1, 10)

    # The core equation is: Total = (Total * Fraction Fallen) + Remaining Pearls
    # Let 'x' be the total number of pearls.
    # x = x * fraction_fallen + remaining_pearls
    # Rearranging to solve for x:
    # x * (1 - fraction_fallen) = remaining_pearls
    # x = remaining_pearls / (1 - fraction_fallen)
    
    total_pearls = int(remaining_pearls / (1 - fraction_fallen))

    print("--- Part 1: How many pearls were there altogether? ---")
    print("Let 'x' be the total number of pearls.")
    print("The equation is formed by summing the parts:")
    # We print the equation with all the numbers as requested.
    print(f"x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + ({eleven} * {eleven} - {seven})")
    print(f"x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + {remaining_pearls}")
    print(f"\nSolving this equation reveals the original number of pearls on the necklace.")
    print(f"Total pearls altogether: {total_pearls}")
    print("\n" + "="*50 + "\n")

    # Part 2: Calculate how many more pearls are needed.

    # Total fallen pearls = Original total - Pearls remaining on the string.
    fallen_pearls = total_pearls - remaining_pearls

    # They find 1/3 of the fallen pearls.
    found_pearls = int(fallen_pearls / 3)

    # The number of pearls they still need is the amount that is still lost.
    needed_pearls = fallen_pearls - found_pearls

    print("--- Part 2: How many more pearls are needed? ---")
    print(f"First, we calculate the number of fallen pearls:")
    print(f"{total_pearls} (total) - {remaining_pearls} (on string) = {fallen_pearls} (fallen)")
    
    print(f"\nThey find 1/3 of the fallen pearls:")
    print(f"{fallen_pearls} / 3 = {found_pearls} (found)")
    
    print(f"\nThe number of pearls they still need is the number that remains lost:")
    print(f"{fallen_pearls} (fallen) - {found_pearls} (found) = {needed_pearls} (needed)")
    print(f"\nThey will need {needed_pearls} more pearls.")

# Run the solver function
solve_pearl_riddle()
<<<304>>>