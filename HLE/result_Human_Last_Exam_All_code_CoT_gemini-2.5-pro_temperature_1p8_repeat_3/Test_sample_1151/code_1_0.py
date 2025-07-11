import fractions

def solve_pearl_riddle():
    """
    This function solves the pearl riddle step-by-step.
    """
    # --- Part 1: How many pearls were there altogether? ---

    print("--- Part 1: How many pearls were there altogether? ---")

    # Step 1: Calculate pearls remaining on the string.
    # "a seven shy of eleven times eleven"
    pearls_on_string = (11 * 11) - 7
    print(f"First, we calculate the number of pearls left on the string: (11 * 11) - 7 = {pearls_on_string} pearls.")

    # Step 2: Define the fractions of fallen pearls.
    f_floor = fractions.Fraction(1, 6)
    f_bed = fractions.Fraction(1, 5)
    f_woman = fractions.Fraction(1, 3)
    f_lover = fractions.Fraction(1, 10)
    
    # Sum the fractions of fallen pearls.
    total_fallen_fraction = f_floor + f_bed + f_woman + f_lover

    # Step 3: Set up and solve the equation for the total number of pearls 'x'.
    # The fraction of pearls remaining on the string is 1 - total_fallen_fraction.
    # This fraction is equal to the 'pearls_on_string' we calculated.
    fraction_on_string = 1 - total_fallen_fraction
    
    # x * fraction_on_string = pearls_on_string  => x = pearls_on_string / fraction_on_string
    total_pearls = int(pearls_on_string / fraction_on_string)

    print("\nLet 'x' be the total number of pearls. The full equation is:")
    print(f"({f_floor.numerator}/{f_floor.denominator})*x + ({f_bed.numerator}/{f_bed.denominator})*x + ({f_woman.numerator}/{f_woman.denominator})*x + ({f_lover.numerator}/{f_lover.denominator})*x + {pearls_on_string} = x")
    print(f"\nThe fractions of fallen pearls sum to {total_fallen_fraction.numerator}/{total_fallen_fraction.denominator}.")
    print(f"This means the {pearls_on_string} pearls on the string are the remaining {fraction_on_string.numerator}/{fraction_on_string.denominator} of the total.")
    print(f"Solving for the total: x = {pearls_on_string} / ({fraction_on_string.numerator}/{fraction_on_string.denominator}) = {total_pearls}")
    print(f"\nSo, there were {total_pearls} pearls altogether on the necklace.")
    
    print("\n" + "="*50 + "\n")

    # --- Part 2: How many more pearls are needed? ---
    
    print("--- Part 2: How many more pearls are they going to need? ---")

    # Step 1: Calculate the number of fallen pearls.
    fallen_pearls = total_pearls - pearls_on_string
    print(f"Total fallen pearls = {total_pearls} (total) - {pearls_on_string} (on string) = {fallen_pearls} pearls.")

    # Step 2: Calculate how many pearls are still needed.
    # They find 1/3, which means 2/3 are still missing.
    still_needed_pearls = int(fallen_pearls * fractions.Fraction(2, 3))
    print(f"They find 1/3 of the fallen pearls, so they are still missing the other 2/3.")
    print(f"Pearls they still need = (2/3) * {fallen_pearls} = {still_needed_pearls}")

    return total_pearls, still_needed_pearls

# Execute the function and get the final answer.
total, needed = solve_pearl_riddle()
final_answer = f"The necklace had {total} pearls in total, and they still need to find {needed} more."

# The final result in the requested format.
# Although a print statement is used here for demonstration,
# this would be the final return value of a program.
print(f"\n<<<{final_answer}>>>")
