# --- Part 1: How many pearls were there altogether? ---

# Calculate the number of pearls remaining on the string from the riddle's text:
# "a seven shy of eleven times eleven"
remaining_on_string = 11 * 11 - 7

# The problem states the total number of pearls (x) is the sum of its parts.
# The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_on_string
# Combining the fractions: 1/6 + 1/5 + 1/3 + 1/10 = 24/30 = 4/5.
# So, the pearls that fell off make up 4/5 of the total.
# This means the remaining pearls make up 1 - 4/5 = 1/5 of the total.
# Therefore, (1/5) * x = remaining_on_string. We can solve for x.

total_pearls = remaining_on_string * 5

# Calculate each part based on the total to display the full equation
floor = total_pearls / 6
bed = total_pearls / 5
woman_saved = total_pearls / 3
lover_caught = total_pearls / 10

print("To find the total number of pearls, we solve the following equation based on the riddle:")
print(f"Total = (Total/6) + (Total/5) + (Total/3) + (Total/10) + {remaining_on_string}")
print("\nSubstituting the final value for Total ({total_pearls}), we can verify the parts:")
print(f"{total_pearls} = {int(floor)} (floor) + {int(bed)} (bed) + {int(woman_saved)} (woman) + {int(lover_caught)} (lover) + {remaining_on_string} (string)")
verification_sum = int(floor) + int(bed) + int(woman_saved) + int(lover_caught) + remaining_on_string
print(f"{total_pearls} = {verification_sum}")
print(f"\nTherefore, there were {total_pearls} pearls altogether.")

print("\n" + "---" * 10 + "\n")

# --- Part 2: How many more pearls are they gonna need? ---

# The number of fallen pearls is the total minus what's left on the string.
fallen_pearls = total_pearls - remaining_on_string

# They manage to find back 1/3rd of the fallen ones.
found_pearls = fallen_pearls / 3

# The number of pearls needed to restore the necklace to its original state.
needed_pearls = fallen_pearls - found_pearls

print("For the second question: How many more pearls are needed?")
print(f"The original necklace had {total_pearls} pearls.")
print(f"Number of pearls that fell off: {total_pearls} - {remaining_on_string} = {fallen_pearls}")
print(f"They find back 1/3 of the fallen pearls: {fallen_pearls} / 3 = {int(found_pearls)}")
print(f"This means they are still missing {fallen_pearls} - {int(found_pearls)} = {int(needed_pearls)} pearls.")
print(f"\nThey are going to need {int(needed_pearls)} more pearls to complete the necklace.")