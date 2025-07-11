import math

def solve_pearl_riddle():
    """
    Solves the two-part pearl necklace word problem.
    """
    
    # --- Part 1: Calculate the total number of pearls ---
    
    # The number of pearls remaining on the string is "seven shy of eleven times eleven".
    remaining_on_string = (11 * 11) - 7
    
    # The fractions of the total pearls that fell.
    fraction_on_floor = 1/6
    fraction_on_bed = 1/5
    fraction_with_woman = 1/3
    fraction_with_lover = 1/10
    
    # The total fraction of pearls that fell is the sum of these parts.
    total_fallen_fraction = fraction_on_floor + fraction_on_bed + fraction_with_woman + fraction_with_lover
    
    # The total number of pearls (x) can be found with the equation:
    # x = (total_fallen_fraction * x) + remaining_on_string
    # x - (total_fallen_fraction * x) = remaining_on_string
    # x * (1 - total_fallen_fraction) = remaining_on_string
    # x = remaining_on_string / (1 - total_fallen_fraction)
    total_pearls = remaining_on_string / (1 - total_fallen_fraction)
    
    # Since the number of pearls must be a whole number, we round it.
    total_pearls = int(round(total_pearls))

    print("--- Part 1: How many pearls were there altogether? ---")
    print(f"The necklace originally had {total_pearls} pearls.\n")
    
    # As requested, output each number in the final equation.
    pearls_on_floor = int(total_pearls * fraction_on_floor)
    pearls_on_bed = int(total_pearls * fraction_on_bed)
    pearls_with_woman = int(total_pearls * fraction_with_woman)
    pearls_with_lover = int(total_pearls * fraction_with_lover)
    
    print("The total is the sum of all the scattered pearls plus those left on the string:")
    print(f"Equation: {pearls_on_floor} (floor) + {pearls_on_bed} (bed) + {pearls_with_woman} (woman) + {pearls_with_lover} (lover) + {remaining_on_string} (string) = {total_pearls} (Total)")

    print("\n" + "="*70 + "\n")

    # --- Part 2: Calculate how many more pearls are needed ---
    
    # Total fallen pearls is the total minus what's left on the string.
    fallen_pearls = total_pearls - remaining_on_string
    
    # They find 1/3 of the fallen pearls.
    found_pearls = int(fallen_pearls * (1/3))
    
    # The number they still need is the amount that is still lost.
    needed_pearls = fallen_pearls - found_pearls
    
    print("--- Part 2: How many more pearls do they need? ---")
    print(f"Total pearls that fell off the string: {fallen_pearls}")
    print(f"Number of pearls they manage to find (1/3 of fallen): {found_pearls}")
    print(f"The number of pearls they still need to find or replace is {needed_pearls}.")
    
    return needed_pearls

# Run the function and get the final answer.
final_answer = solve_pearl_riddle()
print(f"\nFinal Answer: They need {final_answer} more pearls.")

# The final answer in the requested format.
print(f"<<<{final_answer}>>>")
