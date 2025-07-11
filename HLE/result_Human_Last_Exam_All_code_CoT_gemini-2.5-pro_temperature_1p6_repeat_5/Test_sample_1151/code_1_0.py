from fractions import Fraction

# --- Part 1: How many pearls were there altogether? ---

print("--- Question 1: How many pearls were there altogether? ---")

# Step 1: Define the known quantities from the problem.
f_floor = Fraction(1, 6)
f_bed = Fraction(1, 5)
f_woman = Fraction(1, 3)
f_lover = Fraction(1, 10)

# Calculate the number of pearls remaining on the string from the riddle's description.
pearls_on_string_val = (11 * 11) - 7

print(f"\nFirst, we calculate the number of pearls left on the string:")
print(f"Eleven times eleven is 11 * 11 = {11*11}.")
print(f"A 'seven shy' of that is {11*11} - 7 = {pearls_on_string_val} pearls.")

print("\nNext, we set up an equation where 'x' is the total number of pearls.")
print(f"The equation representing the sum of all parts is:")
print(f"x = ({f_floor.numerator}/{f_floor.denominator})x + ({f_bed.numerator}/{f_bed.denominator})x + ({f_woman.numerator}/{f_woman.denominator})x + ({f_lover.numerator}/{f_lover.denominator})x + {pearls_on_string_val}")

# Step 2: Solve the equation for 'x'.
# To solve, we find what fraction of 'x' the known number of pearls represents.
f_fallen_sum = f_floor + f_bed + f_woman + f_lover
f_on_string = 1 - f_fallen_sum

# The number of pearls on the string is equal to this remaining fraction times the total.
# pearls_on_string_val = f_on_string * x
# So, x = pearls_on_string_val / f_on_string
total_pearls = pearls_on_string_val / f_on_string

print("\nTo solve for x, we combine the fractions:")
print(f"Total fraction of fallen pearls = {f_floor} + {f_bed} + {f_woman} + {f_lover} = {f_fallen_sum}")
print(f"The fraction of pearls remaining on the string is 1 - {f_fallen_sum} = {f_on_string}.")
print(f"Therefore, ({f_on_string}) * x = {pearls_on_string_val}.")
print(f"x = {pearls_on_string_val} / {f_on_string} = {int(total_pearls)}")
print(f"\nAnswer 1: There were {int(total_pearls)} pearls on the necklace altogether.")


# --- Part 2: How many more pearls do they need? ---

print("\n\n--- Question 2: How many more pearls are needed? ---")

# Step 3: Calculate the number of pearls needed to complete the necklace again.
fallen_pearls = total_pearls - pearls_on_string_val
found_pearls = fallen_pearls / 3
needed_pearls = fallen_pearls - found_pearls

print(f"\nFirst, we find the number of fallen pearls:")
print(f"{int(total_pearls)} (total) - {pearls_on_string_val} (on string) = {int(fallen_pearls)} fallen pearls.")

print("\nThey manage to find back 1/3 of the fallen ones:")
print(f"(1/3) * {int(fallen_pearls)} = {int(found_pearls)} pearls found.")

print("\nThe number of pearls they still need is the amount that fell minus the amount they found:")
print(f"{int(fallen_pearls)} - {int(found_pearls)} = {int(needed_pearls)}")
print(f"\nAnswer 2: They will need {int(needed_pearls)} more pearls.")