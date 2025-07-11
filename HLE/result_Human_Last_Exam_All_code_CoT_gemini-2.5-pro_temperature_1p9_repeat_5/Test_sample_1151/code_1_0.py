import fractions

# --- Part 1: How many pearls were there altogether? ---
print("--- Part 1: Finding the total number of pearls on the necklace ---")

# Let 'x' be the total number of pearls.
# The number of pearls remaining on the string is "a seven shy of eleven times eleven".
pearls_remained_op_1 = 11
pearls_remained_op_2 = 11
pearls_remained_op_3 = 7
pearls_remained = (pearls_remained_op_1 * pearls_remained_op_2) - pearls_remained_op_3

# The fractions of the total pearls ('x') that fell are:
f_floor = fractions.Fraction(1, 6)   # 1/6 fell to the floor
f_bed = fractions.Fraction(1, 5)     # 1/5 fell upon the bed
f_woman = fractions.Fraction(1, 3)   # 1/3 were saved by the woman
f_lover = fractions.Fraction(1, 10)  # 1/10 were caught by her lover

# The equation is: Total_Pearls = (Sum of Fallen Fractions * Total_Pearls) + Pearls_Remained
# We print the full equation with the numbers from the poem.
print("The equation to solve is:")
print(f"x = ({f_floor.numerator}/{f_floor.denominator})x + ({f_bed.numerator}/{f_bed.denominator})x + ({f_woman.numerator}/{f_woman.denominator})x + ({f_lover.numerator}/{f_lover.denominator})x + (({pearls_remained_op_1} * {pearls_remained_op_2}) - {pearls_remained_op_3})")
print(f"This simplifies to: x = (Sum of Fallen Fractions) * x + {pearls_remained}")

# To solve for x, we rearrange the equation: x - (Sum of Fallen Fractions) * x = 114
# x * (1 - Sum of Fallen Fractions) = 114
# x = 114 / (1 - Sum of Fallen Fractions)

f_fallen_sum = f_floor + f_bed + f_woman + f_lover
f_on_string = 1 - f_fallen_sum
total_pearls = int(pearls_remained / f_on_string)

print(f"\nSolving this equation gives the total number of pearls.")
print(f"Total number of pearls altogether: {total_pearls}")

print("\n" + "="*60 + "\n")

# --- Part 2: How many more pearls are they gonna need? ---
print("--- Part 2: Finding the number of pearls they need to buy ---")

# First, find the total number of pearls that fell off the string.
fallen_pearls = total_pearls - pearls_remained

# They find 1/3rd of the fallen pearls.
fraction_found = fractions.Fraction(1, 3)
found_pearls = int(fallen_pearls * fraction_found)

# The number they need is the difference between what fell and what was found.
needed_pearls = fallen_pearls - found_pearls

print(f"Total pearls that fell: {total_pearls} (total) - {pearls_remained} (remained) = {fallen_pearls} pearls.")
print(f"Pearls they found back (1/3 of fallen): {fraction_found.numerator}/{fraction_found.denominator} * {fallen_pearls} = {found_pearls} pearls.")
print(f"The number of pearls they still need to buy: {fallen_pearls} (fallen) - {found_pearls} (found) = {needed_pearls} pearls.")
