# Set the dimension of the system
d = 3

# Calculate the number of contacts for a fully jammed particle
jammed_contacts = 2 * d

# Calculate the number of contacts for the unstable particle
unstable_contacts = d + 1

# Print the comparison and the name of the unstable particles
print(
    f"In a hard-spheres system with d={d}, "
    f"most particles are jammed with {2} * {d} = {jammed_contacts} contacts.\n"
    f"The unstable particles, which have {d} + {1} = {unstable_contacts} contacts, are known as rattlers."
)