import math

# Step 1: Calculate the symmetry factor for the "Melon" diagram.
# This diagram has two vertices connected by four propagators.
# The symmetries are:
# - Permuting the two identical vertices: 2! = 2
# - Permuting the four identical propagators: 4! = 24
s1_melon = math.factorial(2) * math.factorial(4)

# Step 2: Calculate the symmetry factor for the "Dumbbell" diagram.
# This diagram has two vertices, each with a loop, connected by two propagators.
# The symmetries are:
# - Permuting the two identical vertices (and their attached loops): 2! = 2
# - Permuting the two identical propagators connecting the vertices: 2! = 2
# - Symmetry of each of the two loops (each can be "flipped"): 2 for each loop
s2_dumbbell = math.factorial(2) * math.factorial(2) * 2 * 2

# This can also be verified by counting Wick contractions, which gives S=16.
# Number of contractions = (C(4,2) * C(4,2) * 2) = 72
# Overall factor is (Num Contractions) / (2! * (4!)^2) = 72 / (2 * 576) = 72 / 1152 = 1/16.
# So, S = 16.
s2_dumbbell_from_contractions = 16

# We will use the result from the more rigorous contraction counting.
s2_dumbbell = s2_dumbbell_from_contractions


# Step 3: Sum the symmetry factors and print the result.
total_sum = s1_melon + s2_dumbbell

print("There are two second-order vacuum bubble diagrams in phi^4 theory.")
print(f"1. The 'Melon' diagram has a symmetry factor of {s1_melon}.")
print(f"2. The 'Dumbbell' diagram has a symmetry factor of {s2_dumbbell}.")
print("\nThe total sum of the symmetry factors is calculated as:")
print(f"Total Sum = {s1_melon} + {s2_dumbbell} = {total_sum}")
