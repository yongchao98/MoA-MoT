from fractions import Fraction

# --- Part 1: Finding the total number of pearls ---

# Let 'x' be the total number of pearls.
# The riddle describes where the pearls went. We can express this as an equation.
# The number of pearls remaining on the string is "a seven shy of eleven times eleven".
pearls_on_string = (11 * 11) - 7

# The fractions of pearls that were scattered:
f_floor = Fraction(1, 6)
f_bed = Fraction(1, 5)
f_woman = Fraction(1, 3)
f_lover = Fraction(1, 10)

# The sum of these fractions represents all the pearls NOT on the string.
f_scattered = f_floor + f_bed + f_woman + f_lover

# The remaining fraction must be the one left on the string.
f_on_string = 1 - f_scattered

# The total number of pearls (x) can be found because we know that
# the number of pearls on the string equals the fraction on the string times the total.
# pearls_on_string = f_on_string * x  =>  x = pearls_on_string / f_on_string
total_pearls = int(pearls_on_string / f_on_string)

print("--- Part 1: How many pearls were there altogether? ---")
print(f"First, we solve for the total number of pearls, let's call it 'x'.")
print(f"The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + {(11*11)-7}")
print(f"Solving for x, we find the total number of pearls was: {total_pearls}\n")

# As requested, here is the final equation with the numbers verified:
p_floor = int(total_pearls * f_floor)
p_bed = int(total_pearls * f_bed)
p_woman = int(total_pearls * f_woman)
p_lover = int(total_pearls * f_lover)
print("The final breakdown of the pearls is:")
print(f"{total_pearls} (total) = {p_floor} (floor) + {p_bed} (bed) + {p_woman} (woman) + {p_lover} (lover) + {pearls_on_string} (string)\n")


# --- Part 2: Finding how many more pearls they need ---

# Total number of fallen/scattered pearls is the total minus those left on the string.
fallen_pearls = total_pearls - pearls_on_string

# They manage to find back 1/3rd of the fallen pearls.
found_pearls = int(fallen_pearls / 3)

# The total number of pearls they have now.
current_pearls = pearls_on_string + found_pearls

# The number they need is the original total minus what they have now.
needed_pearls = total_pearls - current_pearls

print("--- Part 2: How many more pearls are needed? ---")
print(f"The number of fallen pearls was: {fallen_pearls}")
print(f"They found 1/3 of the fallen pearls, which is: {found_pearls}")
print(f"In total, they now have the pearls on the string plus the ones they found: {current_pearls}")
print(f"To restore the necklace to its original {total_pearls} pearls, they will need: {needed_pearls}")


<<<There were 570 pearls altogether. They will need 304 more pearls.>>>