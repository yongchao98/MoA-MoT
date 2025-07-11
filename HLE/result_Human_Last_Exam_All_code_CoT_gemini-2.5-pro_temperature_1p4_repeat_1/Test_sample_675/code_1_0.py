# Define the dimension of the system
d = 3

# Calculate the average number of contacts for a particle in the rigid
# backbone of a jammed system (isostatic condition).
z_jammed_val = 2 * d

# Calculate the typical number of contacts for an unstable particle,
# known as a rattler.
z_rattler_val = d + 1

print(f"In a {d}-dimensional hard-sphere system at jamming:")
print("Most particles are stable and have a number of contacts 'z' given by the isostatic condition z = 2d.")
print(f"For d = {d}, the equation is: z = 2 * {d} = {z_jammed_val}")
print("\nHowever, some particles are not strictly jammed and can move locally.")
print("These are called rattlers, and they have a typical number of contacts equal to d+1.")
print(f"For d = {d}, the equation is: z_rattler = {d} + 1 = {z_rattler_val}")