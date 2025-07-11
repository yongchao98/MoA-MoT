# Define the dimension of the system
d = 3

# Calculate the number of contacts for a stable, isostatic particle
# In a jammed system, most particles are isostatic, meaning they have just enough
# contacts to be mechanically stable. For frictionless hard spheres, this is 2*d.
z_stable = 2 * d

# Calculate the number of contacts for the unstable particles described.
# These particles have fewer contacts than required for isostaticity.
z_unstable = d + 1

print(f"In a {d}-dimensional hard sphere system at jamming:")
print(f"The average number of contacts for a stable particle is 2 * {d} = {z_stable}.")
print(f"The number of contacts for the unstable particles is {d} + 1 = {z_unstable}.")
print("\nParticles with fewer contacts than the isostatic limit (z < 2d) are not part of the rigid, stress-bearing backbone of the material.")
print("They are trapped by their neighbors but can still move or 'rattle' within their cages.")
print("\nThe name of these particles is: Rattlers")