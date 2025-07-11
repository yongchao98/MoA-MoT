def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle based on the provided text.
    """

    # Part 1: How many pearls were there altogether?

    # Calculate the number of pearls remaining on the string from the description:
    # "a seven shy of eleven times eleven"
    pearls_on_string = (11 * 11) - 7

    # Let 'X' be the total number of pearls.
    # The problem can be written as an equation where X is the sum of all its parts:
    # X = (1/6)X + (1/5)X + (1/3)X + (1/10)X + pearls_on_string
    #
    # To solve for X, we first isolate X:
    # X - (1/6 + 1/5 + 1/3 + 1/10)X = pearls_on_string
    #
    # Summing the fractions: 1/6 + 1/5 + 1/3 + 1/10 = 5/30 + 6/30 + 10/30 + 3/30 = 24/30 = 4/5
    # So, the equation becomes:
    # X * (1 - 4/5) = pearls_on_string
    # X * (1/5) = pearls_on_string
    # X = pearls_on_string * 5
    
    total_pearls = pearls_on_string * 5

    print("--- Part 1: Total Pearls on the Necklace ---")
    print("The number of pearls left on the string is calculated from the phrase 'a seven shy of eleven times eleven':")
    print(f"(11 * 11) - 7 = 121 - 7 = {pearls_on_string}")
    
    print("\nThe full equation for the total number of pearls (X) is:")
    # Printing each number in the equation as requested
    print(f"X = (1/6)X + (1/5)X + (1/3)X + (1/10)X + {pearls_on_string}")
    
    print("\nBy solving this, we find the total number of pearls.")
    print(f"The total number of pearls on the necklace was: {total_pearls}")

    # Part 2: How many more are they now gonna need?

    # First, find the number of pearls that fell off the string.
    fallen_pearls = total_pearls - pearls_on_string

    # They find back 1/3 of the fallen pearls.
    found_pearls = fallen_pearls / 3

    # The number they still need is the number of fallen pearls they did NOT find (2/3).
    needed_pearls = fallen_pearls - found_pearls
    
    print("\n--- Part 2: Pearls Still Needed ---")
    print(f"The number of fallen pearls is {total_pearls} - {pearls_on_string} = {int(fallen_pearls)}.")
    print(f"They find 1/3rd of the fallen pearls, which is {int(fallen_pearls)} / 3 = {int(found_pearls)}.")
    print(f"The number of pearls still missing is {int(fallen_pearls)} - {int(found_pearls)} = {int(needed_pearls)}.")
    print(f"\nTherefore, the number of additional pearls they need is: {int(needed_pearls)}")

solve_pearl_riddle()
<<<304>>>