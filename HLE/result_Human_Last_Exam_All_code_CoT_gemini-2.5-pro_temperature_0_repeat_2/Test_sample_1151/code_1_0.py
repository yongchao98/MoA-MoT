import fractions

def solve_pearl_riddle():
    """
    This function solves the pearl riddle step-by-step.
    """
    # --- Part 1: How many pearls were there altogether? ---

    # Calculate the number of pearls remaining on the string from the clue:
    # "a seven shy of eleven times eleven"
    pearls_on_string = (11 * 11) - 7

    # The fractions of the total pearls that fell off
    f_floor = fractions.Fraction(1, 6)
    f_bed = fractions.Fraction(1, 5)
    f_woman = fractions.Fraction(1, 3)
    f_lover = fractions.Fraction(1, 10)

    # The sum of these fractions represents the total portion of pearls that fell
    sum_of_fallen_fractions = f_floor + f_bed + f_woman + f_lover

    # The remaining portion is what's left on the string
    fraction_on_string = 1 - sum_of_fallen_fractions

    # The total number of pearls (x) can be found by:
    # pearls_on_string = fraction_on_string * x
    # So, x = pearls_on_string / fraction_on_string
    total_pearls = int(pearls_on_string / fraction_on_string)

    print("--- Part 1: Finding the Total Number of Pearls ---")
    print(f"The number of pearls remaining on the string is (11 * 11) - 7 = {pearls_on_string}.")
    print("The equation for the total number of pearls (x) is:")
    print(f"x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + {pearls_on_string}\n")
    print(f"By solving this, we find the total number of pearls was: {total_pearls}\n")

    # To fulfill the request to show the final equation with numbers:
    # We calculate each part based on the total number of pearls.
    floor_pearls = int(total_pearls * f_floor)
    bed_pearls = int(total_pearls * f_bed)
    woman_pearls = int(total_pearls * f_woman)
    lover_pearls = int(total_pearls * f_lover)
    
    print("Verification of the total with each part calculated:")
    print(f"{total_pearls} = {floor_pearls} (floor) + {bed_pearls} (bed) + {woman_pearls} (woman) + {lover_pearls} (lover) + {pearls_on_string} (string)")
    print(f"Sum of parts: {floor_pearls + bed_pearls + woman_pearls + lover_pearls + pearls_on_string}\n")

    # --- Part 2: How many more pearls are needed? ---

    # The target for the new necklace mentioned in the first line
    target_necklace_size = 500

    # Total pearls that fell off the string
    fallen_pearls = total_pearls - pearls_on_string

    # They find 1/3 of the fallen pearls
    found_pearls = int(fallen_pearls / 3)

    # Total pearls they have now (on string + found)
    total_pearls_possessed = pearls_on_string + found_pearls

    # Pearls they still need to reach the target of 500
    needed_pearls = target_necklace_size - total_pearls_possessed

    print("--- Part 2: Calculating the Number of Pearls Needed ---")
    print(f"Total pearls that fell from the string: {total_pearls} - {pearls_on_string} = {fallen_pearls}.")
    print(f"They find 1/3rd of the fallen pearls: {fallen_pearls} / 3 = {found_pearls}.")
    print(f"Total pearls they now have: {pearls_on_string} (on string) + {found_pearls} (found) = {total_pearls_possessed}.")
    print(f"For a {target_necklace_size} pearl necklace, they still need: {target_necklace_size} - {total_pearls_possessed} = {needed_pearls}.")
    
    return needed_pearls

# Execute the function and store the final answer
final_answer = solve_pearl_riddle()
# The final answer is printed by the function itself.
# The value is returned to be used in the final formatted output.
# print(f"\n<<< {final_answer} >>>") # This line is for the final wrapper, not for the code block itself.

if __name__ == '__main__':
    solve_pearl_riddle()