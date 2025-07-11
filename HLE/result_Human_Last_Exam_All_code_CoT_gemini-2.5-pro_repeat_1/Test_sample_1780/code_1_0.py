import math

# Based on the analysis, the problem is equivalent to finding the number of
# ways to arrange the 3 letters {L, N, S}, since W must always be at the end.
num_items_to_permute = 3

# The number of ways to arrange n distinct items is n! (n factorial).
result = math.factorial(num_items_to_permute)

print(f"The number of valid arrangements is the number of permutations of the {num_items_to_permute} letters L, N, and S.")
print("The formula for this is 3! (3 factorial).")
print(f"The calculation is: 3 * 2 * 1 = {result}")
print("The final number of ways is:")
print(result)