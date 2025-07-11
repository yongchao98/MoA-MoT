# Set the dimension of the hard-sphere system as described in the problem
d = 3

# Calculate the average number of contacts for a particle in the jammed state
# This is given by the isostatic condition, z = 2d
z_jammed_particles = 2 * d

# Calculate the number of contacts for the special, unstable particles
# This is given as d+1 in the problem description
z_unstable_particles = d + 1

print(f"The dimension of the system is d = {d}.")
print(f"The average number of contacts for a typical jammed particle is calculated as: 2 * {d} = {z_jammed_particles}")
print(f"The number of contacts for the special particles described is calculated as: {d} + 1 = {z_unstable_particles}")