# Define the dimension of the system as given in the problem
d = 3

# Define the numerical factors for the contact equations
jammed_factor = 2
unstable_addend = 1

# Calculate the number of contacts for a typical jammed particle (isostatic condition)
# z = 2 * d
z_jammed = jammed_factor * d

# Calculate the number of contacts for the unstable particles
# z = d + 1
z_unstable = d + unstable_addend

# Print the analysis
print(f"In a hard spheres system of dimension d={d}:")
print(f"A strictly jammed particle is expected to have z = {jammed_factor} * {d} = {z_jammed} contacts.")
print(f"The particles in question have z = {d} + {unstable_addend} = {z_unstable} contacts.")
print("\nParticles in a jammed packing that have fewer contacts than required for mechanical stability (i.e., they are not part of the force-bearing structure) and are free to move within their local cage are called 'Rattlers'.")
