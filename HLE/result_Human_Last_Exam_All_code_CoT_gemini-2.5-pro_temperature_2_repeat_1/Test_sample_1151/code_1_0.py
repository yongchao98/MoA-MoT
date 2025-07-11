import math

def solve_pearl_riddle():
    """
    This function solves the pearl necklace riddle based on the provided text.
    """

    # --- Part 1: Find the total number of pearls ---

    print("Step 1: Finding the total number of pearls.")

    # The problem can be represented by the equation:
    # x = x/6 + x/5 + x/3 + x/10 + remaining_pearls
    # where x is the total number of pearls.
    # To solve for x, we rearrange the equation:
    # x - x/6 - x/5 - x/3 - x/10 = remaining_pearls
    # (1 - 1/6 - 1/5 - 1/3 - 1/10) * x = remaining_pearls

    # First, calculate the number of pearls remaining on the string.
    rem_term1 = 11
    rem_term2 = 11
    rem_term3 = 7
    remaining_pearls = (rem_term1 * rem_term2) - rem_term3
    
    print(f"The number of pearls remaining on the string is ({rem_term1} * {rem_term2}) - {rem_term3} = {remaining_pearls}.")

    # The coefficient of x is (1 - 1/6 - 1/5 - 1/3 - 1/10).
    # The common denominator is 30.
    # (30 - 5 - 6 - 10 - 3) / 30 = 6/30 = 1/5.
    # So, (1/5) * x = remaining_pearls
    x_coefficient = 1 - (1/6 + 1/5 + 1/3 + 1/10)
    total_pearls = remaining_pearls / x_coefficient
    
    # Calculate each part for the final equation output
    on_floor = total_pearls / 6
    on_bed = total_pearls / 5
    woman_saved = total_pearls / 3
    lover_caught = total_pearls / 10

    print(f"Solving the equation for the total number of pearls gives: {int(total_pearls)}")
    
    print("\nThe full equation with the solution is:")
    # The prompt requested to output each number in the final equation
    print(
        f"{int(total_pearls)} (Total) = "
        f"{int(on_floor)} (1/6 on floor) + "
        f"{int(on_bed)} (1/5 on bed) + "
        f"{int(woman_saved)} (1/3 saved by woman) + "
        f"{int(lover_caught)} (1/10 caught by lover) + "
        f"{int(remaining_pearls)} (on string)"
    )

    # --- Part 2: Find how many more they need ---
    
    print("\nStep 2: Finding how many more pearls are needed.")
    
    # Calculate the total number of fallen pearls
    fallen_pearls = total_pearls - remaining_pearls
    print(f"Total fallen pearls = {int(total_pearls)} - {int(remaining_pearls)} = {int(fallen_pearls)}")
    
    # Calculate the number of pearls they find back
    found_pearls = fallen_pearls / 3
    print(f"They find back 1/3rd of the fallen pearls: {int(fallen_pearls)} / 3 = {int(found_pearls)}")

    # Calculate the number of pearls still needed
    needed_pearls = fallen_pearls - found_pearls
    print(f"\nThe number of pearls they still need to find is {int(fallen_pearls)} - {int(found_pearls)} = {int(needed_pearls)}.")

    return needed_pearls

# Execute the function and store the final answer
final_answer = solve_pearl_riddle()
print(f"\n<<<304>>>")
