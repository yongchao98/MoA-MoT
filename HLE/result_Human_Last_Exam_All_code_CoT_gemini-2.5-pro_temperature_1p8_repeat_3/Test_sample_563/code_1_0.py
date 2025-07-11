# The number of isomorphism classes of automorphism groups of compact,
# connected Riemann surfaces X for a given genus g.

# For genus g=2, the number of groups is 12.
# This is a known result from the classification of automorphism groups
# of hyperelliptic Riemann surfaces.
num_groups_g2 = 12

# For genus g=3, the number of groups is 36.
# This consists of 35 non-trivial groups and the trivial group
# for a generic surface.
num_groups_g3 = 36

# For genus g=4, the number of groups is 23.
# This consists of 22 non-trivial groups and the trivial group
# for a generic surface.
num_groups_g4 = 23

# The problem asks for the result in the format [g2, g3, g4].
# Let's print the breakdown first, as requested.
print(f"Number of automorphism group classes for genus g=2: {num_groups_g2}")
print(f"Number of automorphism group classes for genus g=3: {num_groups_g3}")
print(f"Number of automorphism group classes for genus g=4: {num_groups_g4}")

# Now, we assemble the final list.
result = [num_groups_g2, num_groups_g3, num_groups_g4]

# Print the final result in the specified format.
print("\nFinal Answer:")
print(result)
