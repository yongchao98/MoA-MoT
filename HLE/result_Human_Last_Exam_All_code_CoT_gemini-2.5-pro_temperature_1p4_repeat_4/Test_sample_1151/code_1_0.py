import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle based on the provided text.
    """

    # Part 1: How many pearls were there altogether?

    # Step 1: Calculate the number of pearls remaining on the string.
    # "a seven shy of eleven times eleven" means (11 * 11) - 7.
    pearls_on_string = 11 * 11 - 7

    # Step 2: Set up and solve the equation for the total number of pearls (x).
    # The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string
    # To solve for x, we rearrange it:
    # x - (1/6 + 1/5 + 1/3 + 1/10)x = pearls_on_string
    # x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = pearls_on_string
    # x = pearls_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))

    sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10
    total_pearls = pearls_on_string / (1 - sum_of_fractions)

    # The number of pearls must be a whole number, so we round it.
    total_pearls = int(round(total_pearls))

    print("--- Part 1: How many pearls were there altogether? ---")

    # To fulfill the request to print each number in the final equation,
    # we calculate each part based on the total.
    pearls_on_floor = int(round(total_pearls * (1/6)))
    pearls_on_bed = int(round(total_pearls * (1/5)))
    pearls_saved_by_woman = int(round(total_pearls * (1/3)))
    pearls_caught_by_lover = int(round(total_pearls * (1/10)))

    # Step 3: Print the complete equation and the answer.
    print("The final equation is the sum of all parts:")
    print(f"{total_pearls} (Total) = {pearls_on_floor} (1/6 fell to the floor) + "
          f"{pearls_on_bed} (1/5 fell on the bed) + "
          f"{pearls_saved_by_woman} (1/3 saved by the woman) + "
          f"{pearls_caught_by_lover} (1/10 caught by the lover) + "
          f"{pearls_on_string} (remained on the string)")
    print(f"\nTherefore, there were {total_pearls} pearls altogether.\n")


    # Part 2: How many more pearls do they need?

    # Step 4: Calculate the total number of fallen pearls.
    fallen_pearls = total_pearls - pearls_on_string

    # Step 5: Calculate how many they find back and how many are still needed.
    found_pearls = math.floor(fallen_pearls / 3)
    needed_pearls = fallen_pearls - found_pearls

    print("--- Part 2: How many more pearls do they need? ---")
    print(f"Total fallen pearls: {fallen_pearls}")
    print(f"Pearls they find back (1/3 of fallen): {found_pearls}")
    print(f"Number of pearls they still need: {needed_pearls}")


solve_pearl_riddle()
<<<304>>>