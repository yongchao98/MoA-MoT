# Define the dimension of the system from the problem description
d = 3

# Calculate the number of contacts for a typical jammed particle
z_jammed = 2 * d

# Calculate the number of contacts for the unstable particle
z_unstable = d + 1

# Print the information in a clear format, showing the calculation
print(f"In d={d}, a typical jammed particle has 2 * {d} = {z_jammed} contacts.")
print(f"The unstable particle described has {d} + 1 = {z_unstable} contacts.")
