import sys
from io import StringIO

def solve_pearl_riddle():
    """
    This function solves the pearl riddle step-by-step.
    """
    # --- Part 1: How many pearls were there altogether? ---
    
    # Step 1: Calculate the number of pearls left on the string
    # "seven shy of eleven times eleven"
    pearls_left_on_string = 11 * 11 - 7

    print("--- Part 1: Finding the total number of pearls ---")
    print("First, we calculate the number of pearls that remained on the string:")
    print(f"11 * 11 - 7 = {pearls_left_on_string}")
    
    # Step 2: Solve for the total number of pearls (x)
    # The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_left_on_string
    # We can rearrange it to: x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = pearls_left_on_string
    # Which simplifies to: x * (1 - 1/6 - 1/5 - 1/3 - 1/10) = pearls_left_on_string
    
    # The sum of the fractions is: 1/6 + 1/5 + 1/3 + 1/10 = 5/30 + 6/30 + 10/30 + 3/30 = 24/30 = 4/5
    # So the equation is: x * (1 - 4/5) = 114  or  x * (1/5) = 114
    
    coefficient = 1 - (1/6 + 1/5 + 1/3 + 1/10)
    total_pearls = pearls_left_on_string / coefficient

    print("\nTo find the total number of pearls (x), we solve the equation:")
    print(f"x = (x/6) + (x/5) + (x/3) + (x/10) + {pearls_left_on_string}")
    print("This simplifies to the final calculation:")
    print(f"x = {pearls_left_on_string} / (1 - 1/6 - 1/5 - 1/3 - 1/10)")
    print(f"Total number of pearls = {int(total_pearls)}")
    
    print("\n----------------------------------------------------\n")

    # --- Part 2: How many more pearls are needed? ---

    print("--- Part 2: Finding how many more pearls are needed ---")
    target_necklace_size = 500
    
    # Step 1: Calculate the number of fallen pearls
    fallen_pearls = total_pearls - pearls_left_on_string
    print(f"1. Total fallen pearls = {int(total_pearls)} (total) - {pearls_left_on_string} (on string) = {int(fallen_pearls)}")
    
    # Step 2: Calculate how many they find back
    found_pearls = fallen_pearls / 3
    print(f"2. They find 1/3 of the fallen pearls: {int(fallen_pearls)} / 3 = {int(found_pearls)}")
    
    # Step 3: Calculate how many pearls they have now
    current_pearls_owned = pearls_left_on_string + found_pearls
    print(f"3. Total pearls they have now = {pearls_left_on_string} (on string) + {int(found_pearls)} (found) = {int(current_pearls_owned)}")
    
    # Step 4: Calculate how many more they need
    pearls_needed = target_necklace_size - current_pearls_owned
    print(f"\n4. To get a {target_necklace_size}-pearl necklace, they need:")
    print(f"{target_necklace_size} (target) - {int(current_pearls_owned)} (owned) = {int(pearls_needed)}")
    print(f"\nSo, they will need {int(pearls_needed)} more pearls for their desired necklace.")

    # Redirect stdout to capture the final answer
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()

    print(int(pearls_needed))

    sys.stdout = old_stdout
    final_answer = captured_output.getvalue().strip()
    return final_answer

final_answer = solve_pearl_riddle()
print(f"\n<<<{final_answer}>>>")