# The dimension of the system
d = 3

# Calculate the number of contacts for a stable, jammed particle
jammed_contacts = 2 * d

# Calculate the number of contacts for the unstable particle
unstable_contacts = d + 1

# The name of these particles
particle_name = "rattlers"

# Print the final result including the calculations
print(f"In a {d}D system, most jammed particles have 2*d = 2*{d} = {jammed_contacts} contacts.")
print(f"The unstable particles, called {particle_name}, have d+1 = {d}+1 = {unstable_contacts} contacts.")