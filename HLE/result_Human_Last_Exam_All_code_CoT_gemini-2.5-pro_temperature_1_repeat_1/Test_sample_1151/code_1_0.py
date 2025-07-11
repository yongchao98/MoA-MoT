import sys
from fractions import Fraction

# This script solves the pearl necklace riddle.

# --- Part 1: How many pearls were there altogether? ---

# Calculate the number of pearls remaining on the string: "a seven shy of eleven times eleven".
remaining_pearls = 11 * 11 - 7

# The fractions of the total pearls that were scattered.
fractions_list = [Fraction(1, 6), Fraction(1, 5), Fraction(1, 3), Fraction(1, 10)]
sum_of_fractions = sum(fractions_list)

# The equation is: Total = (sum_of_fractions * Total) + remaining_pearls
# Rearranging to solve for Total: Total * (1 - sum_of_fractions) = remaining_pearls
# So, Total = remaining_pearls / (1 - sum_of_fractions)
remaining_fraction = 1 - sum_of_fractions
total_pearls = int(remaining_pearls / remaining_fraction)

print("--- Solving the Pearl Riddle ---")
print("\nPart 1: How many pearls were there altogether?")
print(f"First, we calculate the number of pearls remaining on the string: (11 * 11) - 7 = {remaining_pearls}")

print("\nLet 'x' be the total number of pearls. The riddle gives us the equation:")
# The prompt requires printing each number in the final equation.
# The numbers are the coefficients and the constant term.
f1, f2, f3, f4 = "1/6", "1/5", "1/3", "1/10"
print(f"x = ({f1})x + ({f2})x + ({f3})x + ({f4})x + {remaining_pearls}")

print(f"\nBy solving this equation, we find the total number of pearls was: {total_pearls}")

# --- Part 2: How many more pearls are needed? ---

print("\n----------------------------------------")
print("\nPart 2: How many more pearls are needed?")

# Calculate the total number of fallen pearls.
fallen_pearls = total_pearls - remaining_pearls
print(f"Total number of fallen pearls: {total_pearls} (total) - {remaining_pearls} (on string) = {fallen_pearls}")

# They find 1/3 of the fallen pearls.
found_pearls = fallen_pearls / 3
print(f"They find 1/3 of the fallen pearls: {fallen_pearls} / 3 = {int(found_pearls)}")

# Calculate the total number of pearls they have now.
current_pearls_count = remaining_pearls + found_pearls

# Calculate how many more pearls they need to complete the necklace.
needed_pearls = total_pearls - current_pearls_count
print(f"Pearls they have now: {remaining_pearls} (on string) + {int(found_pearls)} (found) = {int(current_pearls_count)}")
print(f"\nTo complete the necklace, they still need: {int(needed_pearls)} pearls.")

# Final answer in the required format
# Suppress the final answer from the regular output to avoid duplication
original_stdout = sys.stdout
sys.stdout = open('/dev/null', 'w')
print(f'<<<{int(needed_pearls)}>>>')
sys.stdout.close()
sys.stdout = original_stdout
