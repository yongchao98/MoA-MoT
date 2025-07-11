# --- Solving the Pearl Riddle ---

# Part 1: How many pearls were there altogether?
print("Part 1: How many pearls were there altogether?")

# Step 1: Calculate the number of pearls remaining on the string.
# "a seven shy of eleven times eleven"
remained_on_string = (11 * 11) - 7
print(f"First, we calculate the pearls remaining on the string: (11 * 11) - 7 = {remained_on_string}")

# Step 2: Solve for the total number of pearls.
# Let 'x' be the total number of pearls. The equation is:
# x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remained_on_string
# Rearranging the equation to solve for x:
# x - (1/6 + 1/5 + 1/3 + 1/10)x = remained_on_string
# The sum of fractions is 24/30 or 4/5. So, (1 - 4/5)x = remained_on_string
# (1/5)x = remained_on_string
# x = remained_on_string * 5
total_pearls = remained_on_string * 5
print(f"Solving the equation for the total number of pearls, we find the answer is {int(total_pearls)}.")

# Step 3: Verify the result by showing the full equation with the numbers.
floor_pearls = total_pearls * (1/6)
bed_pearls = total_pearls * (1/5)
woman_pearls = total_pearls * (1/3)
lover_pearls = total_pearls * (1/10)
print("\nVerification of the equation:")
print(f"Total Pearls = (Fell to floor) + (Fell on bed) + (Woman saved) + (Lover caught) + (Remained on string)")
print(f"{int(total_pearls)} = {int(floor_pearls)} + {int(bed_pearls)} + {int(woman_pearls)} + {int(lover_pearls)} + {int(remained_on_string)}")
print(f"The sum is indeed {int(floor_pearls + bed_pearls + woman_pearls + lover_pearls + remained_on_string)}.")
print(f"\nTherefore, there were {int(total_pearls)} pearls altogether.")


# Part 2: How many more pearls are needed?
print("\n------------------------------------")
print("\nPart 2: How many more pearls are they gonna need?")

# Step 1: Calculate the total number of fallen pearls.
total_fallen = total_pearls - remained_on_string
print(f"The total number of fallen pearls is: {int(total_pearls)} - {int(remained_on_string)} = {int(total_fallen)}")

# Step 2: Calculate how many they find back (1/3 of fallen).
found_back = total_fallen * (1/3)
print(f"They manage to find back 1/3rd of the fallen ones: {int(total_fallen)} / 3 = {int(found_back)}")

# Step 3: Calculate how many are still missing and therefore needed.
needed_more = total_fallen - found_back
print(f"The number of pearls still missing is: {int(total_fallen)} - {int(found_back)} = {int(needed_more)}")
print(f"\nTherefore, they will need {int(needed_more)} more pearls.")
