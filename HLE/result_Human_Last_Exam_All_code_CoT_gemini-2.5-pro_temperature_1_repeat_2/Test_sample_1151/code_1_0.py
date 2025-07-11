import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle based on the provided text.
    """
    # Part 1: How many pearls were there altogether?

    # Calculate the number of pearls remaining on the string
    # "seven shy of eleven times eleven"
    pearls_on_string = (11 * 11) - 7

    # The equation is T = T/6 + T/5 + T/3 + T/10 + pearls_on_string
    # T - (T/6 + T/5 + T/3 + T/10) = pearls_on_string
    # T * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = pearls_on_string
    # T = pearls_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))
    
    sum_of_fractions = (1/6) + (1/5) + (1/3) + (1/10)
    total_pearls = pearls_on_string / (1 - sum_of_fractions)
    
    # Ensure the result is an integer, as we can't have partial pearls
    total_pearls = int(round(total_pearls))

    # Calculate the value of each part for the final equation printout
    p1 = int(total_pearls / 6)
    p2 = int(total_pearls / 5)
    p3 = int(total_pearls / 3)
    p4 = int(total_pearls / 10)

    print("--- Part 1: Total Pearls ---")
    print(f"First, we calculate the pearls remaining on the string: (11 * 11) - 7 = {pearls_on_string}")
    print("We solve the equation: Total = (Total/6) + (Total/5) + (Total/3) + (Total/10) + Pearls_on_string")
    print("\nThe final equation with the calculated numbers is:")
    print(f"{total_pearls} = {p1} (fell to floor) + {p2} (fell on bed) + {p3} (woman saved) + {p4} (lover caught) + {pearls_on_string} (on string)")
    print(f"\nTherefore, there were {total_pearls} pearls altogether.")

    # Part 2: How many more pearls are needed?
    
    # "fallen ones" are all pearls not on the string
    fallen_pearls = total_pearls - pearls_on_string
    
    # They find back 1/3rd of the fallen ones
    found_pearls = math.floor(fallen_pearls / 3)
    
    # Total they have now = what was on the string + what they found
    current_pearls = pearls_on_string + found_pearls
    
    # Needed pearls = original total - what they have now
    needed_pearls = total_pearls - current_pearls
    
    print("\n--- Part 2: Needed Pearls ---")
    print(f"The total number of fallen pearls is {total_pearls} - {pearls_on_string} = {fallen_pearls}.")
    print(f"They find back 1/3 of these, which is {fallen_pearls} / 3 = {found_pearls}.")
    print(f"The number of pearls they need to restore the necklace is {total_pearls} - ({pearls_on_string} + {found_pearls}) = {needed_pearls}.")
    print(f"\nSo, they will need {needed_pearls} more pearls.")

solve_pearl_riddle()
<<<570>>>