import math

# This problem can be solved by finding the number of valid permutations
# of the sequence of reflections.

# The multiset of mirrors for the reflections is {G1, G1, G2, G3, G4}.
# Total number of reflections is 5.
num_reflections_g1 = 2
num_reflections_other = 3
total_reflections = num_reflections_g1 + num_reflections_other

print(f"The sequence of reflections involves {num_reflections_g1} reflections on G1 and 1 on each of the other {num_reflections_other} mirrors.")
print(f"The total number of reflections in the path is {total_reflections}.")
print("-" * 30)

# Step 1: Calculate the total number of permutations of the multiset {G1, G1, G2, G3, G4}.
# The formula is n! / k!, where n is the total number of items and k is the number of repetitions.
total_perms_numerator = math.factorial(total_reflections)
total_perms_denominator = math.factorial(num_reflections_g1)
total_permutations = total_perms_numerator // total_perms_denominator

print("Step 1: Calculate the total number of permutations if there were no restrictions.")
print(f"Total permutations = {total_reflections}! / {num_reflections_g1}! = {total_perms_numerator} / {total_perms_denominator} = {total_permutations}")
print("-" * 30)

# Step 2: Calculate the number of invalid permutations where the two G1s are adjacent.
# We treat the two G1s as a single block: {'G1G1', G2, G3, G4}.
# The number of items to permute is now (total_reflections - 1).
num_items_for_invalid_perm = total_reflections - 1
invalid_permutations = math.factorial(num_items_for_invalid_perm)

print("Step 2: Calculate the number of invalid permutations (where G1 reflections are adjacent).")
print(f"Treating the two G1 reflections as a single block, we find the permutations of {num_items_for_invalid_perm} items.")
print(f"Invalid permutations = {num_items_for_invalid_perm}! = {invalid_permutations}")
print("-" * 30)

# Step 3: Subtract the invalid from the total to get the number of valid ways.
valid_ways = total_permutations - invalid_permutations

print("Step 3: Calculate the final number of ways by subtracting invalid from total.")
print(f"Number of ways = {total_permutations} - {invalid_permutations} = {valid_ways}")

print("<<<36>>>")