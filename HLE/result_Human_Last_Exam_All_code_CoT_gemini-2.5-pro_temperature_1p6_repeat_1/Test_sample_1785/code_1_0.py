# Step 1: Define the minimum length of the components for the simplest
# nontrivial 3-component link.

# The simplest nontrivial 3-component link consists of a minimal 2-component link (a Hopf link)
# and a minimal unlinked component (an unknot).

# According to results in lattice knot theory, the minimal length of a Hopf link is 12 edges,
# realized as two interlocking loops of 6 edges each.
len_c1 = 6  # Length of the first component of the Hopf link
len_c2 = 6  # Length of the second component of the Hopf link

# The minimal length for any single loop (an unknot) on a lattice is a 1x1 square.
len_c3 = 4  # Length of the third, unlinked component

# Step 2: Calculate the total minimum number of edges by summing the lengths of the three components.
total_length = len_c1 + len_c2 + len_c3

# Step 3: Print the final equation and the result.
# The final equation shows the contribution of each component to the total length.
print(f"{len_c1} + {len_c2} + {len_c3} = {total_length}")