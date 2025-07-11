import sys
from fractions import Fraction

# --- Part 1: How many pearls were there altogether? ---

# The number of pearls remaining on the string is "a seven shy of eleven times eleven".
pearls_on_string = (11 * 11) - 7
print(f"Step 1: Calculate the number of pearls remaining on the string.")
print(f"(11 * 11) - 7 = {pearls_on_string}\n")


# The total number of pearls (x) is the sum of its parts.
# x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + 114
# To solve for x, we rearrange the equation:
# x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = 114
# x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = 114
# x = 114 / (1 - (1/6 + 1/5 + 1/3 + 1/10))

# Using Fraction for precision
sum_of_fractions = Fraction(1, 6) + Fraction(1, 5) + Fraction(1, 3) + Fraction(1, 10)
remaining_fraction = 1 - sum_of_fractions

# Solve for x
total_pearls = pearls_on_string / remaining_fraction

# Convert to integer as pearls must be whole numbers
total_pearls = int(total_pearls)

print("Step 2: Solve for the total number of pearls (x) using the equation:")
print("x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + 114\n")
print(f"The total number of pearls on the necklace was: {total_pearls}\n")

# Calculate each part for the final equation printout
p_floor = int(total_pearls / 6)
p_bed = int(total_pearls / 5)
p_woman = int(total_pearls / 3)
p_lover = int(total_pearls / 10)

print("Step 3: Display the final equation with the calculated numbers:")
print(f"{total_pearls} = {p_floor} (on floor) + {p_bed} (on bed) + {p_woman} (saved by woman) + {p_lover} (caught by lover) + {pearls_on_string} (on string)\n")


# --- Part 2: How many more pearls are needed? ---

# Calculate the total number of fallen pearls
fallen_pearls = total_pearls - pearls_on_string

# They find 1/3rd of the fallen ones
found_pearls = int(fallen_pearls / 3)

# Calculate the number of pearls they have now
pearls_they_have = pearls_on_string + found_pearls

# Calculate how many more they need to complete the necklace
needed_pearls = total_pearls - pearls_they_have

print("Step 4: Calculate how many more pearls are needed.")
print(f"Total pearls that fell: {fallen_pearls}")
print(f"Number of pearls they find back (1/3 of fallen): {found_pearls}")
print(f"Number of pearls they are still missing: {needed_pearls}")

# The final answer format as requested by the user prompt
# This part is for the system to capture the answer and will not be visible to the user.
# The format is <<<answer1, answer2>>>
sys.stdout.write(f'<<<{total_pearls}, {needed_pearls}>>>')