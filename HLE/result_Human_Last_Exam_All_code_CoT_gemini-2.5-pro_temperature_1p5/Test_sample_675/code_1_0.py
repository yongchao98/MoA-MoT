# Define the dimension of the system
d = 3

# Calculate the average number of contacts for a jammed particle
z_jammed = 2 * d

# Calculate the number of contacts for the unstable particle
z_unstable = d + 1

# Explain the situation with the calculated values and provide the name
print(
    f"In a hard spheres system with dimension d={d}, most jammed particles have on average z = 2*d = {z_jammed} contacts.\n"
    f"The particles in question, which are not strictly jammed and have z = d+1 = {z_unstable} contacts, are called 'rattlers'."
)