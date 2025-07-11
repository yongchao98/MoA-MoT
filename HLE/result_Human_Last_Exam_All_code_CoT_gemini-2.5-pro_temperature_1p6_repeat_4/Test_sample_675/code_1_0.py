# Set the dimension of the system
d = 3

# Calculate the average number of contacts for a jammed particle (isostatic condition)
z_jammed = 2 * d

# Calculate the number of contacts for the unstable, non-jammed particles
z_unstable = d + 1

# Print the calculations and results
print("In a hard sphere system at jamming in d=3:")
print(f"The number of contacts for a typical, jammed particle is z = 2 * d = {2} * {d} = {z_jammed}.")
print(f"The number of contacts for an unstable particle is z = d + 1 = {d} + {1} = {z_unstable}.")
print("\nThe particles with d+1 contacts are not part of the rigid, force-bearing backbone of the structure and can move within their local cage.")
print("These particles are known as 'rattlers'.")
