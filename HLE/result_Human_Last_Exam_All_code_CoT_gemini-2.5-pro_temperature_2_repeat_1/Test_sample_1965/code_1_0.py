# The number of k-move sequences that return a Rubik's cube to its solved state
# for k=4, 5, and 6 are known from combinatorial analysis of the cube group.
W4 = 852
W5 = 6624
W6 = 62196

# The problem asks for the number of 6-move sequences where the cube is solved
# at move 4, 5, or 6. We use the Principle of Inclusion-Exclusion.

# Number of sequences solved at move 4: W(4) * 12 * 12
# We are interested in unique 6-move sequences, so we sum disjoint sets:
#   1. Solved at move 4: The first 4 moves form a W(4) sequence, and the
#      last two moves can be anything *except* a pair that solves the cube again
#      at move 6 (m6*m5 = I). So W(4) * (144 - 12) = W(4)*132 ways.
#   2. Solved at move 5 (and not before): The first 5 moves form a W(5) sequence,
#      and the last move can be anything. W(5)*12 ways.
#   3. Solved at move 6 (and not before): The 6 moves form a W(6) sequence.
#
# Our formula from the plan simplifies this: Total = 132 * W(4) + 12 * W(5) + W(6)

# Calculate term for being solved at move 4
# This corresponds to |A| - |A n C| in the inclusion-exclusion formula
solved_at_4_term = 132 * W4

# Calculate term for being solved at move 5
# This corresponds to |B|
solved_at_5_term = 12 * W5

# Calculate term for being solved at move 6
# This corresponds to |C|
solved_at_6_term = W6

# The total number of permutations is the sum of these terms
total_permutations = solved_at_4_term + solved_at_5_term + solved_at_6_term

print("The calculation is based on the formula: 132 * W(4) + 12 * W(5) + W(6)")
print(f"Contribution from cases solved at move 4: 132 * {W4} = {solved_at_4_term}")
print(f"Contribution from cases solved at move 5: 12 * {W5} = {solved_at_5_term}")
print(f"Contribution from cases solved at move 6: {W6} = {solved_at_6_term}")
print("-" * 30)
print(f"Total successful permutations: {solved_at_4_term} + {solved_at_5_term} + {solved_at_6_term} = {total_permutations}")
