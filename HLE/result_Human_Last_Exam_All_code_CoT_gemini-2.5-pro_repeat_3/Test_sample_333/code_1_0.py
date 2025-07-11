# Define the known variables based on the problem description.
# I is the horizontal distance from the gun to the highest point of elevation.
I = 500

# x1 is the landing distance of the first fragment from the gun.
# To find the maximum distance for the second fragment, we must assume
# the first fragment lands at the minimum possible distance, i.e., at the gun.
x1 = 0

# The center of mass (CM) of the system continues on the original trajectory.
# The landing spot of the CM would be at 2 * I.
# For two equal fragments, the landing spot of the CM is also (x1 + x2) / 2.
# So, 2 * I = (x1 + x2) / 2, which simplifies to 4 * I = x1 + x2.
# From this, we can find the landing spot of the second fragment, x2.
# x2 = 4 * I - x1

# Now we print the final equation with the numbers substituted in.
print(f"The equation relating the fragments' landing spots (x1, x2) and the distance to the apex (I) is: 4 * I = x1 + x2")
print(f"To find the maximum distance, we use I = {I} m and x1 = {x1} m.")
print(f"Substituting these values: 4 * {I} = {x1} + x2")

# Calculate the result
max_distance = 4 * I - x1

print(f"Solving for x2: {4*I} = {x1} + x2")
print(f"The maximum distance from the gun you can safely be is x2 = {max_distance} m.")