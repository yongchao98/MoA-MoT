import fractions

# --- Part 1: How many pearls were there altogether? ---

# 1. Calculate the number of pearls remaining on the string.
remaining_on_string = (11 * 11) - 7

# 2. Define and sum the fractions of pearls that fell.
f_floor = fractions.Fraction(1, 6)
f_bed = fractions.Fraction(1, 5)
f_woman = fractions.Fraction(1, 3)
f_lover = fractions.Fraction(1, 10)
fallen_fraction = f_floor + f_bed + f_woman + f_lover

# 3. Calculate the fraction of pearls that remained on the string.
remaining_fraction = 1 - fallen_fraction

# 4. Calculate the total number of pearls.
# (remaining_fraction) * total_pearls = remaining_on_string
total_pearls = remaining_on_string / remaining_fraction

print("--- Part 1: How many pearls were there altogether? ---")
print(f"The number of pearls remaining on the string is calculated from 'a seven shy of eleven times eleven':")
print(f"  (11 * 11) - 7 = {remaining_on_string}")
print("\nThe full equation for the total number of pearls ('x') is:")
print(f"  x = ({f_floor.numerator}/{f_floor.denominator})x + ({f_bed.numerator}/{f_bed.denominator})x + ({f_woman.numerator}/{f_woman.denominator})x + ({f_lover.numerator}/{f_lover.denominator})x + {remaining_on_string}")
print(f"\nThe pearls on the string represent the fraction of the total that did not fall: 1 - {fallen_fraction.numerator}/{fallen_fraction.denominator} = {remaining_fraction.numerator}/{remaining_fraction.denominator}")
print(f"So, ({remaining_fraction.numerator}/{remaining_fraction.denominator}) of the total pearls equals {remaining_on_string}.")
print(f"Total Pearls = {remaining_on_string} / ({remaining_fraction.numerator}/{remaining_fraction.denominator}) = {int(total_pearls)}")
print(f"\nTherefore, there were {int(total_pearls)} pearls on the necklace altogether.")


# --- Part 2: How many more pearls are needed? ---

# 1. Calculate the number of fallen pearls.
fallen_pearls = total_pearls - remaining_on_string

# 2. Calculate how many fallen pearls are found.
found_pearls = fallen_pearls / 3

# 3. Calculate how many more pearls are needed.
current_pearls_count = remaining_on_string + found_pearls
needed_pearls = total_pearls - current_pearls_count

print("\n\n--- Part 2: How many more pearls are needed? ---")
print(f"The number of pearls that fell was: {int(total_pearls)} - {remaining_on_string} = {int(fallen_pearls)}")
print(f"They manage to find 1/3 of the fallen pearls: {int(fallen_pearls)} / 3 = {int(found_pearls)}")
print(f"Total pearls they have now = {remaining_on_string} (on string) + {int(found_pearls)} (found) = {int(current_pearls_count)}")
print(f"To restore the necklace, they need: {int(total_pearls)} (total) - {int(current_pearls_count)} (current) = {int(needed_pearls)}")
print(f"\nThey will need {int(needed_pearls)} more pearls.")

<<<304>>>