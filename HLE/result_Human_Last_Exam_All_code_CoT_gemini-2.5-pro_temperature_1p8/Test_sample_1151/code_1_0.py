import math

def solve_necklace_riddle():
    """
    Solves the pearl necklace riddle by first determining the total number of pearls
    and then calculating how many more are needed after finding a portion of the lost ones.
    """
    # Part 1: How many pearls were there altogether?

    # Calculate the number of pearls remaining on the string from the riddle's clue.
    # "a seven shy of eleven times eleven"
    remaining_val_1 = 11
    remaining_val_2 = 11
    shy_of_val = 7
    remaining_on_string = (remaining_val_1 * remaining_val_2) - shy_of_val

    # The fractions of the total pearls that fell off.
    f_floor = 1/6
    f_bed = 1/5
    f_woman = 1/3
    f_lover = 1/10

    # The total fraction of pearls that fell off the string.
    total_fallen_fraction = f_floor + f_bed + f_woman + f_lover

    # The fraction of pearls that remained on the string is 1 minus the fallen fraction.
    remaining_fraction = 1 - total_fallen_fraction

    # Calculate the total number of pearls.
    # total_pearls * remaining_fraction = remaining_on_string
    # So, total_pearls = remaining_on_string / remaining_fraction
    total_pearls = int(round(remaining_on_string / remaining_fraction))

    print("--- Part 1: How many pearls were there altogether? ---")
    print(f"The number of pearls remaining on the string is ({remaining_val_1} * {remaining_val_2}) - {shy_of_val} = {remaining_on_string}.")
    print("The equation to solve for the Total Pearls (x) is:")
    print(f"x = (x * {f_floor:.2f}) + (x * {f_bed:.2f}) + (x * {f_woman:.2f}) + (x * {f_lover:.2f}) + {remaining_on_string}")
    print("\nBy solving for x, we find the total number of pearls was {}.\n".format(total_pearls))

    # Print the final equation with all numbers filled in.
    pearls_on_floor = int(total_pearls * f_floor)
    pearls_on_bed = int(total_pearls * f_bed)
    pearls_saved_by_woman = int(total_pearls * f_woman)
    pearls_caught_by_lover = int(total_pearls * f_lover)
    
    print("The final equation with the numbers is:")
    print(f"{total_pearls} = {pearls_on_floor} + {pearls_on_bed} + {pearls_saved_by_woman} + {pearls_caught_by_lover} + {remaining_on_string}")
    print(f"Checking the sum of the parts: {pearls_on_floor + pearls_on_bed + pearls_saved_by_woman + pearls_caught_by_lover + remaining_on_string}")
    

    # Part 2: How many more are they gonna need?
    print("\n--- Part 2: How many more pearls are needed? ---")

    # The number of pearls that fell is the total minus what's left on the string.
    fallen_pearls = total_pearls - remaining_on_string
    print(f"Total number of fallen pearls: {fallen_pearls}")

    # They find back 1/3rd of the fallen ones.
    found_pearls = int(fallen_pearls / 3)
    print(f"Number of pearls they find back (1/3 of fallen): {found_pearls}")

    # The total number of pearls they have now.
    current_pearls_count = remaining_on_string + found_pearls
    print(f"Total pearls they have now (on string + found): {current_pearls_count}")

    # The number of pearls needed to restore the necklace.
    needed_pearls = total_pearls - current_pearls_count
    print(f"\nTo restore the necklace from {current_pearls_count} back to {total_pearls} pearls, they need:")
    print(f"{total_pearls} - {current_pearls_count} = {needed_pearls} more pearls.")
    
    return needed_pearls

# Run the function and capture the final answer.
final_answer = solve_necklace_riddle()