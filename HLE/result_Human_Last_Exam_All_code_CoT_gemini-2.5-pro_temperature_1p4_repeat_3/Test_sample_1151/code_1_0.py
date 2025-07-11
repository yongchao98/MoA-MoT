import math

def solve_pearl_riddle():
    """
    Solves a multi-part word puzzle about a broken pearl necklace.
    """

    # --- Part 1: How many pearls were there altogether? ---

    print("Part 1: How many pearls were there altogether?")
    
    # Calculate the number of pearls remaining on the string
    pearls_on_string = 11 * 11 - 7
    print(f"First, we calculate the pearls remaining on the string from the clue 'a seven shy of eleven times eleven':")
    print(f"11 * 11 - 7 = 121 - 7 = {pearls_on_string} pearls.")

    # The equation for the total number of pearls 'x' is:
    # x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string
    # This simplifies to: x - (1/6 + 1/5 + 1/3 + 1/10)x = pearls_on_string
    # The sum of fractions is 4/5. So: x - (4/5)x = pearls_on_string
    # Which means: (1/5)x = pearls_on_string
    # Therefore: x = pearls_on_string * 5

    total_pearls = pearls_on_string * 5
    
    print(f"\nBy solving the riddle's equation, we find the total number of pearls on the necklace was {total_pearls}.")
    
    # Calculate each fractional part to display the full equation
    floor_pearls = total_pearls / 6
    bed_pearls = total_pearls / 5
    woman_pearls = total_pearls / 3
    lover_pearls = total_pearls / 10
    
    print("\nThe full equation with the calculated total is:")
    # Using math.trunc to ensure we display whole numbers, as you can't have a fraction of a pearl.
    print(f"{total_pearls} (Total) = {math.trunc(floor_pearls)} (Fell to floor) + {math.trunc(bed_pearls)} (Fell on bed) + {math.trunc(woman_pearls)} (Saved by woman) + {math.trunc(lover_pearls)} (Caught by lover) + {pearls_on_string} (Remained on string)")
    
    print("-" * 50)

    # --- Part 2: How many more pearls do they need? ---
    
    print("\nPart 2: How many more pearls are they now gonna need?")
    
    total_fallen = total_pearls - pearls_on_string
    print(f"The total number of fallen pearls is {total_pearls} - {pearls_on_string} = {math.trunc(total_fallen)}.")
    
    found_pearls = total_fallen / 3
    print(f"They manage to find back 1/3rd of the fallen ones: {math.trunc(total_fallen)} / 3 = {math.trunc(found_pearls)} pearls.")
    
    needed_pearls = total_fallen - found_pearls
    print(f"\nTherefore, the number of pearls they still need to find is {math.trunc(total_fallen)} - {math.trunc(found_pearls)} = {math.trunc(needed_pearls)}.")
    
    return math.trunc(needed_pearls)

# Execute the function and get the final answer
final_answer = solve_pearl_riddle()
print(f"\n<<<>>>")
print(f"<<<{final_answer}>>>")