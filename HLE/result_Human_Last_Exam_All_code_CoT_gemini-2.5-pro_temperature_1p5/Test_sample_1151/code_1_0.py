import math

def solve_pearl_riddle():
    """
    Solves the pearl necklace riddle based on the provided story.
    """
    # --- Part 1: Find the total number of pearls ---

    # The number of pearls remaining on the string is "a seven shy of eleven times eleven"
    remaining_on_string = (11 * 11) - 7

    # Let 'x' be the total number of pearls.
    # The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_on_string
    # To solve for x, we rearrange to: x - (1/6 + 1/5 + 1/3 + 1/10)x = remaining_on_string
    # This simplifies to: x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = remaining_on_string
    
    # Sum of the fractional parts that fell off
    sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10
    
    # Solve for x (total_pearls)
    total_pearls_float = remaining_on_string / (1 - sum_of_fractions)
    
    # Since pearls must be a whole number, we round the result.
    total_pearls = int(round(total_pearls_float))

    print("--- Part 1: How many pearls were there altogether? ---")
    print(f"First, we calculate the pearls left on the string: (11 * 11) - 7 = {remaining_on_string}")
    print("\nThe equation to find the total number of pearls (x) is derived from the riddle:")
    # Printing each number in the equation as requested
    print(f"x = (1/6)*x + (1/5)*x + (1/3)*x + (1/10)*x + {remaining_on_string}")
    print(f"\nBy solving this equation, we find the total number of pearls.")
    print(f"Answer: There were {total_pearls} pearls altogether on the necklace.\n")

    # --- Part 2: Find how many more pearls they need ---

    # Calculate the number of pearls that fell off the string
    fallen_pearls = total_pearls - remaining_on_string
    
    # They find 1/3 of the fallen pearls
    found_pearls = int(round(fallen_pearls / 3))
    
    # Calculate how many pearls they have now
    current_pearls = remaining_on_string + found_pearls
    
    # Calculate how many more they need to get back to the original total
    needed_pearls = total_pearls - current_pearls

    print("--- Part 2: How many more pearls are they going to need? ---")
    print(f"The number of pearls that fell off was: {total_pearls} - {remaining_on_string} = {fallen_pearls}")
    print(f"They managed to find 1/3 of the fallen pearls: {fallen_pearls} / 3 = {found_pearls}")
    print(f"So, they now have a total of: {remaining_on_string} (on string) + {found_pearls} (found) = {current_pearls} pearls.")
    print(f"To restore the necklace to its original {total_pearls} pearls, they will need:")
    print(f"{total_pearls} (original) - {current_pearls} (current) = {needed_pearls} more pearls.")

# Execute the function to print the solution
solve_pearl_riddle()