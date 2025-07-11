# Define the dimension of the system as provided in the problem
d = 3

# Calculate the number of contacts for a typical, isostatic particle in the jammed state.
# The formula given is z = 2d.
jammed_contacts = 2 * d

# Calculate the number of contacts for the special, unstable particles.
# The formula given is z = d + 1.
unstable_contacts = d + 1

# Print the calculated values and the context.
# The prompt requires outputting each number in the final equation.
print("In a hard spheres system at jamming:")
print(f"The number of contacts for a fully jammed particle (isostatic condition) is calculated as: 2 * {d} = {jammed_contacts}")
print(f"The number of contacts for the unstable, rattling particle is calculated as: {d} + 1 = {unstable_contacts}")

# Print the name for these specific particles.
print("\nThe name for particles that have fewer contacts than required for isostaticity (d+1 contacts) and can move back-and-forth within their cage is:")
print("rattlers")