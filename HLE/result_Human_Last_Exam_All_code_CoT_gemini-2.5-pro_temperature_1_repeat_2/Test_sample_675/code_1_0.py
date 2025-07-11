# Define the dimension of the system as provided in the problem.
d = 3

# Calculate the average number of contacts for a typically jammed particle.
# The equation is z = 2d.
jammed_contacts = 2 * d

# Calculate the number of contacts for the specific unstable particles.
# The equation is d + 1.
unstable_contacts = d + 1

# Print the information along with the name of the particles.
# The code outputs each number involved in the descriptive equations.
print(f"In a {d}-dimensional system, most jammed particles have 2 * {d} = {jammed_contacts} contacts.")
print(f"The particles in question, which are unstable and can move, have only {d} + 1 = {unstable_contacts} contacts.")
print("These particles are called: rattlers.")