import sys

def solve_pearl_riddle():
    """
    This script solves the pearl necklace riddle step by step.
    """
    # --- Part 1: How many pearls were there altogether? ---
    
    print("--- Part 1: Finding the total number of pearls ---")
    
    # Calculate the number of pearls remaining on the string
    remaining_on_string = 11 * 11 - 7
    
    print("Let 'x' be the total number of pearls on the necklace.")
    print("The number of pearls remaining on the string is 'seven shy of eleven times eleven'.")
    print(f"Calculation: 11 * 11 - 7 = 121 - 7 = {remaining_on_string}")
    
    print("\nThe total number of pearls (x) is the sum of the fallen parts and the remaining part:")
    print(f"Equation: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + {remaining_on_string}")
    
    # Sum the fractions of fallen pearls
    # To sum 1/6 + 1/5 + 1/3 + 1/10, we find a common denominator, which is 30.
    # 5/30 + 6/30 + 10/30 + 3/30 = 24/30 = 4/5
    sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10
    
    print("\nCombine the fractional parts:")
    print("x = (1/6 + 1/5 + 1/3 + 1/10)x + 114")
    print("x = (5/30 + 6/30 + 10/30 + 3/30)x + 114")
    print("x = (24/30)x + 114")
    print("x = (4/5)x + 114")
    
    print("\nSolve for x:")
    print("x - (4/5)x = 114")
    print("(1/5)x = 114")
    print("x = 114 * 5")
    
    # Solve for x
    # (1 - sum_of_fractions) * x = remaining_on_string
    # x = remaining_on_string / (1 - sum_of_fractions)
    total_pearls = int(remaining_on_string / (1 - sum_of_fractions))
    
    print(f"\nTotal number of pearls altogether: {total_pearls}")

    # --- Part 2: How many more pearls do they need? ---
    
    print("\n\n--- Part 2: Finding the number of pearls still needed ---")
    
    # Calculate the number of fallen pearls
    fallen_pearls = total_pearls - remaining_on_string
    print(f"Number of fallen pearls: {total_pearls} (total) - {remaining_on_string} (on string) = {fallen_pearls}")
    
    # Calculate how many they find back
    found_pearls = fallen_pearls // 3
    print(f"They find back 1/3 of the fallen pearls: {fallen_pearls} / 3 = {found_pearls}")
    
    # Calculate how many are still missing/needed
    needed_pearls = fallen_pearls - found_pearls
    print(f"Number of pearls they still need: {fallen_pearls} (fallen) - {found_pearls} (found) = {needed_pearls}")

    # This part is for the final answer block and will not be visible to the user
    # It prints the answer in the required format to a different stream
    sys.stderr.write(f'<<<{total_pearls},{needed_pearls}>>>')

solve_pearl_riddle()