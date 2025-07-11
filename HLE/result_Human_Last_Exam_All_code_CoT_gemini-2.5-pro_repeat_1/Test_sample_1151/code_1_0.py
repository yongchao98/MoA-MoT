import math

def solve_pearl_riddle():
    """
    Solves the pearl riddle by first finding the total number of pearls
    and then calculating how many more are needed.
    """
    # Part 1: How many pearls were there altogether?

    print("--- Part 1: Finding the Total Number of Pearls ---")

    # The number of pearls remaining on the string is "seven shy of eleven times eleven"
    pearls_on_string = 11 * 11 - 7

    # The fractions of the total pearls ('x') that fell off the string
    # x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string
    # To solve for x, we rearrange the equation:
    # x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = pearls_on_string
    # x * (1 - 1/6 - 1/5 - 1/3 - 1/10) = pearls_on_string
    # x = pearls_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))
    sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10
    total_pearls = pearls_on_string / (1 - sum_of_fractions)

    # Since pearls must be whole numbers, we round the result.
    total_pearls = int(round(total_pearls))

    print(f"The number of pearls remaining on the string is (11 * 11) - 7 = {pearls_on_string}.")
    print("The total number of pearls is the sum of the parts that fell and the parts that remained.")

    # Calculate the number of pearls for each fraction for the final equation printout
    p_floor = int(round(total_pearls / 6))
    p_bed = int(round(total_pearls / 5))
    p_woman = int(round(total_pearls / 3))
    p_lover = int(round(total_pearls / 10))

    print("\nThe final equation with the numbers is:")
    print(f"{p_floor} (fell to floor) + {p_bed} (fell on bed) + {p_woman} (saved by woman) + {p_lover} (caught by lover) + {pearls_on_string} (on string) = {total_pearls}")
    print(f"\nThere were {total_pearls} pearls altogether on the necklace.")


    # Part 2: How many more are they now gonna need?

    print("\n--- Part 2: Calculating the Number of Needed Pearls ---")

    # Total fallen pearls is the total minus what's left on the string
    fallen_pearls = total_pearls - pearls_on_string
    print(f"Total fallen pearls: {total_pearls} - {pearls_on_string} = {fallen_pearls}")

    # They find 1/3rd of the fallen ones
    found_pearls = int(round(fallen_pearls / 3))
    print(f"Pearls they find back: {fallen_pearls} / 3 = {found_pearls}")

    # The number they need is the number of fallen pearls they did NOT find
    needed_pearls = fallen_pearls - found_pearls
    print(f"Pearls they still need to find or replace: {fallen_pearls} - {found_pearls} = {needed_pearls}")

    print(f"\nThey will need {needed_pearls} more pearls.")
    
    return needed_pearls

# Execute the function and capture the final answer
final_answer = solve_pearl_riddle()
# The final answer is the result of the second question.
# To follow the requested format, we print it at the very end.
# print(f"<<<{final_answer}>>>") # This is for the final answer format instruction.

if __name__ == '__main__':
    solve_pearl_riddle()