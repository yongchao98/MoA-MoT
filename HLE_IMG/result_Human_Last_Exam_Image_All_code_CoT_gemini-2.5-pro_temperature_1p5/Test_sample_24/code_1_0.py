import sys

# Step 1: Identify the genus of the Heegaard diagram.
# The number of alpha curves (or beta curves) is 3.
g = 3

# Step 2: From the diagram, we observe that the alpha and beta curves are disjoint.
# This implies the relators in the fundamental group presentation are trivial.

# Step 3: Determine the fundamental group.
# For a genus-g Heegaard diagram with trivial relators, the fundamental group
# is the free group on g generators, F_g.
# For g=3, this is F_3.

# Step 4: Identify the corresponding 3-manifold.
# The prime 3-manifold with fundamental group F_g is the connected sum (#)
# of g copies of the product of a 1-sphere and a 2-sphere (S^1 x S^2).

# Step 5: Construct and print the formula for the manifold.
# We explicitly define the numbers used in the equation as requested.
sphere_dim_1 = 1
sphere_dim_2 = 2

# Create the string for one component of the connected sum
component = f"(S^{sphere_dim_1} \u00d7 S^{sphere_dim_2})" # using unicode for 'x'

# Create the full manifold description
# It is the connected sum of 'g' such components.
manifold_formula = f" # ".join([component] * g)

# Print the final answer
print("The Heegaard diagram represents the three-manifold M given by the following equation:")
# The final equation prints each number: the '1' and '2' for the sphere dimensions, repeated 3 times.
print(f"M = {manifold_formula}")
