import sys

# The dimension of the system as given in the problem description
d = 3

# For most particles in a jammed state, the average number of contacts is z = 2d
# This is known as the isostatic condition.
z_jammed = 2 * d

# For the specific, unstable particles mentioned, the number of contacts is d+1
z_unstable = d + 1

# The name for these particles is "rattlers".
particle_name = "rattlers"

# Print the final answer, including the numbers from the problem's equations.
print(f"In a system of dimension d={d}, isostatically jammed particles have, on average, z = 2*d = {z_jammed} contacts.")
print(f"The particles with only d+1 = {d}+1 = {z_unstable} contacts, which are caged but not strictly jammed, are called '{particle_name}'.")
