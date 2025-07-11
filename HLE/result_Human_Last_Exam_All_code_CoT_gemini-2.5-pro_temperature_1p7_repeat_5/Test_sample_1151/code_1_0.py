import math

# Part 1: Find the total number of pearls

# Let X be the total number of pearls. The equation is:
# X = (1/6)X + (1/5)X + (1/3)X + (1/10)X + Remaining
# We first need to calculate the value of 'Remaining'

remaining_on_string = (11 * 11) - 7

# Now, we solve for X.
# X - (1/6)X - (1/5)X - (1/3)X - (1/10)X = remaining_on_string
# To combine the fractions, we find a common denominator, which is 30.
# (30/30)X - (5/30)X - (6/30)X - (10/30)X - (3/30)X = remaining_on_string
# ((30 - 5 - 6 - 10 - 3)/30)X = remaining_on_string
# (6/30)X = remaining_on_string
# (1/5)X = remaining_on_string
# X = remaining_on_string * 5

total_pearls = remaining_on_string * 5

# Calculate the number of pearls for each part of the equation
p_floor = total_pearls * (1/6)
p_bed = total_pearls * (1/5)
p_woman = total_pearls * (1/3)
p_lover = total_pearls * (1/10)

print("The solved equation for the total number of pearls is:")
print(f"{total_pearls} = {int(p_floor)} (floor) + {int(p_bed)} (bed) + {int(p_woman)} (woman) + {int(p_lover)} (lover) + {remaining_on_string} (string)")
print("-" * 20)
print(f"Total number of pearls on the necklace: {total_pearls}")
print("-" * 20)


# Part 2: Calculate how many more pearls they need

# Calculate the total number of fallen pearls
fallen_pearls = total_pearls - remaining_on_string

# Calculate the number of pearls they find back (1/3rd of fallen ones)
found_pearls = fallen_pearls * (1/3)

# Calculate the number of pearls they still need
needed_pearls = fallen_pearls - found_pearls

print(f"Total fallen pearls: {int(fallen_pearls)}")
print(f"Pearls they find back (1/3 of fallen): {int(found_pearls)}")
print(f"Number of pearls they still need: {int(needed_pearls)}")

<<<304>>>