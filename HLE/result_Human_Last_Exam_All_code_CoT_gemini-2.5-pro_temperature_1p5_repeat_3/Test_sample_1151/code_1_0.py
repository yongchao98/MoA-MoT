import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle by setting up and solving an equation,
    then calculates the number of pearls still needed.
    """
    
    # Part 1: How many pearls were there altogether?
    
    # Calculate the number of pearls that remained on the string.
    # "seven shy of eleven times eleven"
    pearls_on_string = (11 * 11) - 7
    
    # The fractions of the total pearls that fell are 1/6, 1/5, 1/3, and 1/10.
    # The equation is: X = (1/6)X + (1/5)X + (1/3)X + (1/10)X + pearls_on_string
    # To solve for X, we rearrange the equation:
    # X - (1/6 + 1/5 + 1/3 + 1/10)X = pearls_on_string
    # X * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = pearls_on_string
    # X = pearls_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))
    
    sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10
    
    # Calculate total_pearls using the rearranged equation
    total_pearls = pearls_on_string / (1 - sum_of_fractions)
    
    # Since pearls are whole numbers, we round and convert to an integer.
    total_pearls = int(round(total_pearls))

    # Calculate the number for each part to display in the final equation.
    pearls_floor = int(round(total_pearls * (1/6)))
    pearls_bed = int(round(total_pearls * (1/5)))
    pearls_woman = int(round(total_pearls * (1/3)))
    pearls_lover = int(round(total_pearls * (1/10)))
    
    print("--- Part 1: Total Pearls on the Necklace ---")
    print("The riddle describes the distribution of pearls from a broken necklace.")
    print("We can find the total number of pearls by setting up an equation where the total (X) equals the sum of its parts.")
    print("\nTotal Pearls = (Fell to floor) + (Fell on bed) + (Saved by woman) + (Caught by lover) + (Remained on string)")
    
    # Print the equation with all the numbers, as requested.
    print("\nThe solved equation with the actual numbers is:")
    print(f"{total_pearls} = {pearls_floor} + {pearls_bed} + {pearls_woman} + {pearls_lover} + {pearls_on_string}")
    
    print(f"\nTherefore, there were {total_pearls} pearls on the necklace altogether.")
    
    print("\n" + "="*50 + "\n")
    
    # Part 2: How many more pearls are needed?
    
    # Total number of pearls that fell off the string.
    fallen_pearls = total_pearls - pearls_on_string
    
    # They find 1/3rd of the fallen ones.
    found_pearls = int(round(fallen_pearls / 3))
    
    # Total pearls they have now is what was on the string plus what they found.
    current_pearls_total = pearls_on_string + found_pearls
    
    # The number needed is the original total minus what they have now.
    needed_pearls = total_pearls - current_pearls_total
    
    print("--- Part 2: How Many More Pearls Are Needed ---")
    print(f"Of the {fallen_pearls} pearls that fell, they managed to find 1/3rd, which is {found_pearls} pearls.")
    print(f"To restore the necklace to its original {total_pearls} pearls, they will need {needed_pearls} more pearls.")

# Run the function to print the solution.
solve_pearl_riddle()