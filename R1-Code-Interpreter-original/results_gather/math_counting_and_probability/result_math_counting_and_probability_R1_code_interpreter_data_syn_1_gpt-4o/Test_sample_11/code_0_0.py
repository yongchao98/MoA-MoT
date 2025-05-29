import math

# Total permutations of 7 children
total_permutations = math.factorial(7)

# Permutations where no two boys are adjacent
daughters_permutations = math.factorial(3)
boys_permutations = math.factorial(4)
no_adjacent_boys = daughters_permutations * boys_permutations

# Calculate the number of ways with at least 2 boys together
at_least_two_boys_together = total_permutations - no_adjacent_boys

print(at_least_two_boys_together)