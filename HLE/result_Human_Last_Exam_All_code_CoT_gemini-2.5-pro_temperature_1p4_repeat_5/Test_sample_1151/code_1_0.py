from fractions import Fraction

# --- Part 1: How many pearls were there altogether? ---

# First, calculate the number of pearls that remained on the string from the clue:
# "a seven shy of eleven times eleven"
pearls_on_string = 11 * 11 - 7

# Define the fractions of the total pearls that fell off the string
frac_floor = Fraction(1, 6)
frac_bed = Fraction(1, 5)
frac_woman = Fraction(1, 3)
frac_lover = Fraction(1, 10)

# The total fraction of pearls that fell is the sum of these parts
sum_of_fallen_fractions = frac_floor + frac_bed + frac_woman + frac_lover

# The fraction of pearls that remained on the string is 1 minus the fallen fraction
fraction_on_string = 1 - sum_of_fallen_fractions

# The total number of pearls can be found because we know the number of pearls
# that corresponds to the 'fraction_on_string'
# total_pearls * fraction_on_string = pearls_on_string
total_pearls = pearls_on_string / fraction_on_string

print("--- Part 1: How many pearls were there altogether? ---")
print("Let 'x' be the total number of pearls on the necklace.")

print("\nFirst, we calculate the number of pearls that remained on the string:")
print(f"Calculation for 'a seven shy of eleven times eleven': 11 * 11 - 7 = {pearls_on_string}")

print("\nNext, we set up the equation based on the parts of the necklace:")
print(f"x = ({frac_floor})x + ({frac_bed})x + ({frac_woman})x + ({frac_lover})x + {pearls_on_string}")

print("\nTo solve for x, we group the 'x' terms:")
print(f"x - ({frac_floor} + {frac_bed} + {frac_woman} + {frac_lover})x = {pearls_on_string}")
print(f"x * (1 - {sum_of_fallen_fractions}) = {pearls_on_string}")
print(f"x * ({fraction_on_string}) = {pearls_on_string}")
print(f"x = {pearls_on_string} / {fraction_on_string}")
print(f"x = {int(total_pearls)}")

print(f"\nTherefore, there were {int(total_pearls)} pearls altogether on the necklace.")


# --- Part 2: How many more pearls are they gonna need? ---

# Calculate the number of pearls that fell off
fallen_pearls = total_pearls - pearls_on_string

# Calculate how many of the fallen pearls they find
found_pearls = fallen_pearls / 3

# Calculate the final number of pearls they still need
needed_pearls = total_pearls - (pearls_on_string + found_pearls)

print("\n--- Part 2: How many more pearls are they gonna need? ---")
print(f"\nThe total number of pearls that fell is the original total minus those still on the string:")
print(f"Fallen pearls = {int(total_pearls)} - {pearls_on_string} = {int(fallen_pearls)}")

print("\nThey manage to find back 1/3rd of the fallen ones:")
print(f"Found pearls = {int(fallen_pearls)} / 3 = {int(found_pearls)}")

print("\nThe number of pearls they still need is the original total minus all the pearls they have now:")
print(f"Needed pearls = {int(total_pearls)} - ({pearls_on_string} + {int(found_pearls)})")
print(f"Needed pearls = {int(total_pearls)} - {int(pearls_on_string + found_pearls)}")
print(f"Needed pearls = {int(needed_pearls)}")

print(f"\nTherefore, they are going to need {int(needed_pearls)} more pearls.")