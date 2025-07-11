import math

# Step 1: Calculate the total number of distributions (S)
# This is the number of ways to arrange 25 items consisting of 5 groups of 5 identical items.
# S = 25! / (5! * 5! * 5! * 5! * 5!)
N = 25
num_types = 5
items_per_type = 5

# Numerator for S is N!
s_numerator = math.factorial(N)

# Denominator for S is (5!)^5
s_denominator = math.factorial(items_per_type) ** num_types

S = s_numerator // s_denominator

# Step 2: Calculate the number of favorable distributions (F)
# A favorable distribution requires each person to have a unique "dominant" type.
# This is only possible if each person has all 5 items of their assigned dominant type.
# F is the number of ways to assign the 5 types to the 5 people, which is 5!.
F = math.factorial(num_types)

# Step 3: Calculate the probability P = F / S
# The probability is F/S. The numbers F and S are integers.
# The result will be presented as an equation.

print("Step 1: The total number of ways to distribute the items is S.")
print(f"S = 25! / (5!^5)")
print(f"S = {s_numerator} / ({math.factorial(items_per_type)}^{num_types})")
print(f"S = {S}\n")

print("Step 2: The number of favorable distributions is F.")
print(f"F = 5!")
print(f"F = {F}\n")

print("Step 3: The probability P is F / S.")
print(f"P = F / S = {F} / {S}")

# The problem is a theoretical one, so the final answer is the fraction itself.
# No floating point calculation is needed to preserve precision.
# The final result is the equation F/S with the computed integer values.
# Final output should be in the format P = F / S
print("\nThe final equation is:")
print(f"P = {F} / {S}")
