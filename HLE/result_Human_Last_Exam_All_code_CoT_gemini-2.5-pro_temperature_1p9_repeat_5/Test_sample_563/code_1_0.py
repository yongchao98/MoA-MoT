# The task is to find the number of isomorphism classes of automorphism groups
# for compact, connected Riemann surfaces of genus g = 2, 3, and 4.
# These values are not computed from a simple formula but are deep results from
# the classification theory of Riemann surfaces and their automorphism groups.
# We will provide the known values based on this mathematical literature.

# For a Riemann surface of genus g=2, it is known that there are 12
# non-isomorphic groups that can arise as the full automorphism group.
num_groups_g2 = 12

# For a Riemann surface of genus g=3, the classification is more involved,
# leading to a larger set of 36 possible non-isomorphic automorphism groups.
num_groups_g3 = 36

# For a Riemann surface of genus g=4, the number of non-isomorphic
# automorphism groups is 23.
num_groups_g4 = 23

# We will now print these numbers in the requested list format.
# The f-string explicitly includes each number in the final printed output.
print(f"[{num_groups_g2},{num_groups_g3},{num_groups_g4}]")
