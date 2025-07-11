import math

# Given parameters
N = 25  # Total number of items
T = 5   # Number of item types (and individuals)
C = 5   # Copies of each item type

# Step 1: Calculate S, the total number of ways to distribute the items.
# This is the multinomial coefficient C(25; 5, 5, 5, 5, 5), which represents
# the total number of distinct permutations of the 25 items.
S_numerator = math.factorial(N)
S_denominator = math.factorial(C) ** T
S = S_numerator // S_denominator

# Step 2: Calculate F, the number of favorable distributions.
# A favorable outcome is where each individual has a unique item type for which
# they have strictly more items than anyone else.
# We consider the case where each individual j receives all 5 items of a unique type T_j.
# The number of ways to assign the 5 unique types to the 5 individuals is 5!.
# Each specific assignment corresponds to exactly one arrangement of the 25 items.
# So, the number of favorable outcomes is 5!.
F = math.factorial(T)

# Step 3: Output the probability P = F / S.
# The problem asks to output the numbers in the final equation.
print("Given:")
print(f"- N = {N} items")
print(f"- T = {T} types of items, with {C} copies of each type")
print(f"- Each of the {T} individuals receives {C} items.")
print("\nCalculation:")
print(f"S = {N}! / ({C}!)**{T} = {S}")
print(f"F = {T}! = {F}")
print("\nThe probability P = F / S is:")
print(f"P = {F} / {S}")