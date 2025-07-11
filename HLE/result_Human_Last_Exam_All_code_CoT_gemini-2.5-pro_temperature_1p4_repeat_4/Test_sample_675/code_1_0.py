import sys

# The dimension of the system
d = 3

# Calculate the average number of contacts for a stable, jammed particle
# z = 2d
z_jammed_factor_1 = 2
z_jammed = z_jammed_factor_1 * d

# Calculate the number of contacts for the unstable particles in question
# z = d+1
z_unstable_factor_1 = 1
z_unstable = d + z_unstable_factor_1

print(f"In a hard spheres system with dimension d = {d}:")
print(f"Most jammed particles have, on average, z = {z_jammed_factor_1} * {d} = {z_jammed} contacts.")
print(f"Particles that are not strictly jammed and have z = {d} + {z_unstable_factor_1} = {z_unstable} contacts are called 'rattlers'.")
