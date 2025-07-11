import math

# --- Part 1: Find the total number of pearls ---

# Step 1: Calculate the number of pearls remaining on the string
remained_on_string = 11 * 11 - 7

# Step 2: Solve for the total number of pearls (x)
# The equation is: x = x/6 + x/5 + x/3 + x/10 + remained_on_string
# This can be rearranged to: x - (x/6 + x/5 + x/3 + x/10) = remained_on_string
# The sum of the fractions is (1/6 + 1/5 + 1/3 + 1/10) = 24/30 = 4/5
# So, x - (4/5)x = remained_on_string
# which simplifies to (1/5)x = remained_on_string
# Therefore, x = remained_on_string * 5
total_pearls = remained_on_string * 5

# Calculate the number of pearls for each fraction for the final equation printout
on_floor = total_pearls / 6
on_bed = total_pearls / 5
saved_by_woman = total_pearls / 3
caught_by_lover = total_pearls / 10

print("Step 1: Finding the total number of pearls on the necklace.")
print(f"The number of pearls that remained on the string is: 11 * 11 - 7 = {remained_on_string}")
print("The equation to find the total number of pearls (x) is: x = x/6 + x/5 + x/3 + x/10 + pearls_remaining")
print("Substituting the final values into the equation:")
print(f"{int(total_pearls)} = {int(on_floor)} (1/6) + {int(on_bed)} (1/5) + {int(saved_by_woman)} (1/3) + {int(caught_by_lover)} (1/10) + {remained_on_string}")
print(f"So, the total number of pearls was {int(total_pearls)}.\n")

# --- Part 2: Find how many more pearls are needed ---

# Step 3: Calculate the total number of fallen pearls
fallen_pearls = total_pearls - remained_on_string

# Step 4: Calculate the number of pearls found and still needed
found_pearls = fallen_pearls / 3
needed_pearls = fallen_pearls - found_pearls

print("Step 2: Finding how many more pearls are needed.")
print(f"Total number of pearls that fell from the necklace: {int(total_pearls)} - {remained_on_string} = {int(fallen_pearls)}")
print(f"Number of pearls they find back (1/3 of fallen): {int(fallen_pearls)} / 3 = {int(found_pearls)}")
print(f"Number of pearls they still need: {int(fallen_pearls)} - {int(found_pearls)} = {int(needed_pearls)}")

# The final answer is the number of pearls they still need
final_answer = int(needed_pearls)
print(f"\n<<< {final_answer} >>>")