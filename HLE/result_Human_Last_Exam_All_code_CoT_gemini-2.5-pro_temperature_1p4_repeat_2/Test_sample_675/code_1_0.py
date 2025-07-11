import math

# The dimension of the system
d = 3

# Calculate the average number of contacts for a jammed particle
# The formula is z = 2 * d
z_jammed = 2 * d

# Calculate the number of contacts for the unstable particles
# The formula is d + 1
z_unstable = d + 1

print(f"The system is in d = {d} dimensions.")
print(f"The average number of contacts for most jammed particles is z = 2 * {d} = {z_jammed}.")
print(f"The number of contacts for the special, unstable particles is {d} + 1 = {z_unstable}.")
print("\nThese unstable particles, which are not strictly jammed and have a characteristic rattling motion, are called:")
print("rattlers")
