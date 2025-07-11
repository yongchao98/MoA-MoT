def solve_pearl_riddle():
    """
    This function solves the pearl riddle step-by-step.
    """
    # --- Part 1: How many pearls were there altogether? ---
    
    print("--- Part 1: Calculating the total number of pearls ---")
    
    # "a seven shy of eleven times eleven pearls remained upon the string"
    remaining_on_string_val1 = 11
    remaining_on_string_val2 = 11
    remaining_on_string_val3 = 7
    remaining_on_string = (remaining_on_string_val1 * remaining_on_string_val2) - remaining_on_string_val3
    
    print("Let 'x' be the total number of pearls on the necklace.")
    print("The fallen pearls are fractions of the total: 1/6, 1/5, 1/3, and 1/10.")
    print(f"The number of pearls remaining on the string is ({remaining_on_string_val1} * {remaining_on_string_val2}) - {remaining_on_string_val3} = {remaining_on_string}.")
    print("\nSo, the full equation is:")
    # We must output each number in the equation, so we print it out explicitly
    frac_1 = 1
    frac_2 = 6
    frac_3 = 1
    frac_4 = 5
    frac_5 = 1
    frac_6 = 3
    frac_7 = 1
    frac_8 = 10
    print(f"x = ({frac_1}/{frac_2})x + ({frac_3}/{frac_4})x + ({frac_5}/{frac_6})x + ({frac_7}/{frac_8})x + {remaining_on_string}")

    print("\nTo solve for x, we rearrange the equation:")
    print(f"x - (({frac_1}/{frac_2}) + ({frac_3}/{frac_4}) + ({frac_5}/{frac_6}) + ({frac_7}/{frac_8}))x = {remaining_on_string}")
    print("Combining the fractions (1/6 + 1/5 + 1/3 + 1/10) gives 24/30, which simplifies to 4/5.")
    print(f"x - (4/5)x = {remaining_on_string}")
    print(f"(1/5)x = {remaining_on_string}")
    
    # Solve for x: x = remaining_on_string * 5
    solve_val_1 = 1
    solve_val_2 = 5
    total_pearls = remaining_on_string * solve_val_2

    print(f"x = {remaining_on_string} * {solve_val_2} / {solve_val_1}")
    print(f"Total number of pearls altogether: {total_pearls}\n")

    # --- Part 2: How many more do they need? ---

    print("--- Part 2: Calculating the number of pearls they still need ---")

    # Fallen pearls = total - remaining
    fallen_pearls = total_pearls - remaining_on_string
    print(f"First, we find the number of pearls that fell: {total_pearls} (total) - {remaining_on_string} (on string) = {fallen_pearls} fallen pearls.")
    
    # They find 1/3rd of the fallen pearls
    find_frac_numerator = 1
    find_frac_denominator = 3
    found_pearls = fallen_pearls / find_frac_denominator
    # Ensure it's an integer for clear output
    found_pearls = int(found_pearls)
    
    print(f"They find {find_frac_numerator}/{find_frac_denominator} of the fallen pearls: {fallen_pearls} / {find_frac_denominator} = {found_pearls} pearls found.")
    
    # Needed pearls = fallen - found
    needed_pearls = fallen_pearls - found_pearls
    print(f"The number of pearls they still need is: {fallen_pearls} (fallen) - {found_pearls} (found) = {needed_pearls}.")
    
    return total_pearls, needed_pearls

# Run the function and capture the final answer
total, needed = solve_pearl_riddle()
final_answer = f"Total pearls: {total}, Pearls still needed: {needed}"

# Per instructions, print the final answer in the specified format
print(f"\n<<<Total pearls: {total}, Pearls still needed: {needed}>>>")
