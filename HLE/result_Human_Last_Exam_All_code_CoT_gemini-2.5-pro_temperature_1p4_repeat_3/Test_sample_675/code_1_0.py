# Define the dimension of the system from the problem description
d = 3

# Calculate the number of contacts for a typically jammed, stable particle
# according to the isostatic condition z = 2d
num_contacts_jammed = 2 * d

# Calculate the number of contacts for the special, unstable particles
# as described in the problem
num_contacts_unstable = d + 1

# Print the context and the final answer
print(f"In a system of dimension d={d}:")
print(f"A typically jammed particle is isostatic and has z = 2 * d = 2 * {d} = {num_contacts_jammed} contacts.")
print(f"The special particles described are under-constrained and have z = d + 1 = {d} + 1 = {num_contacts_unstable} contacts.")
print("\nThese particles, which are not strictly jammed and exhibit a back-and-forth movement, are called:")
print("Rattlers")