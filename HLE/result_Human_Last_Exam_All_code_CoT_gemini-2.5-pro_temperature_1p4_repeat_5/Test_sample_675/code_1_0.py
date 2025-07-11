# Define the dimension of the system
d = 3

# Calculate the average number of contacts for a jammed particle
# This is based on Maxwell's criterion for mechanical stability, z = 2d
z_jammed = 2 * d

# Calculate the number of contacts for the unstable particles as described
z_unstable = d + 1

# Print the explanation and the answer
print(
    f"In a hard-sphere system in dimension d={d}, most jammed particles have, on average, "
    f"z = 2*d = {z_jammed} contacts.\n"
    f"The particles that are not strictly jammed and have an unstable, "
    f"weak contact with only d+1 = {z_unstable} neighbors are called 'rattlers'."
)
