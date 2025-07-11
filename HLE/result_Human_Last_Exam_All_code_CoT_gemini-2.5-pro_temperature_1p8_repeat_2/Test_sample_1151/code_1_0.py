import math

# --- Part 1: How many pearls were there altogether? ---

# First, calculate the number of pearls that remained on the string.
# The riddle describes this as "a seven shy of eleven times eleven".
remaining_pearls = (11 * 11) - 7

# Let 'x' be the total number of pearls. The equation based on the riddle is:
# x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_pearls
# To solve for x, we can rearrange the equation:
# x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = remaining_pearls
# Factor out x:
# x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = remaining_pearls
# The sum of the fractions is 24/30, or 4/5.
# x * (1 - 4/5) = remaining_pearls
# x * (1/5) = remaining_pearls
# x = remaining_pearls * 5

total_pearls = remaining_pearls * 5

print("--- Question 1: How many pearls were there altogether? ---")
print("The number of pearls remaining on the string is calculated as:")
print(f"Remaining Pearls = (11 * 11) - 7 = {remaining_pearls}")
print("\nThe final equation to find the total number of pearls (x) is:")
print(f"x = {remaining_pearls} * 5")
print(f"Total Pearls = {int(total_pearls)}")
print(f"\nThere were {int(total_pearls)} pearls on the necklace altogether.\n")


# --- Part 2: How many more pearls are they now gonna need? ---

# First, find the total number of pearls that had fallen off.
fallen_pearls = total_pearls - remaining_pearls

# They manage to find back 1/3rd of the fallen ones.
found_pearls = fallen_pearls / 3

# The total number of pearls they currently have.
current_pearls_total = remaining_pearls + found_pearls

# The number of pearls they still need to get back to the original total.
needed_pearls = total_pearls - current_pearls_total


print("--- Question 2: How many more pearls do they need? ---")
print("To find how many more pearls they need, we follow these steps:")
print(f"1. Pearls that fell off: {int(total_pearls)} (Total) - {remaining_pearls} (Remaining) = {int(fallen_pearls)}")
print(f"2. Pearls they found (1/3 of fallen): {int(fallen_pearls)} / 3 = {int(found_pearls)}")
print(f"3. Total pearls they have now: {remaining_pearls} (Remaining) + {int(found_pearls)} (Found) = {int(current_pearls_total)}")
print("4. To restore the necklace, the number of pearls they still need is:")
print(f"Needed Pearls = {int(total_pearls)} (Total) - {int(current_pearls_total)} (Have Now) = {int(needed_pearls)}")

print(f"\nThey will need {int(needed_pearls)} more pearls.")