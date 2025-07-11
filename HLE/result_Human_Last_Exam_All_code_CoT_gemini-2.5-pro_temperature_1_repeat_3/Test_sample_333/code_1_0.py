# Define the known variables based on the problem description.
# I is the horizontal distance from the gun to the highest point of elevation.
I = 500  # in meters

# According to the principle of conservation of momentum, the center of mass
# of the fragments will continue on the original trajectory of the projectile.
# The total range of the original projectile would have been 2 * I.
# This is the landing position of the center of mass (R_cm).
R_cm = 2 * I

# The projectile splits into two equal parts.
# The position of the center of mass is the average of the positions of the two fragments.
# R_cm = (R1 + R2) / 2
# where R1 and R2 are the landing distances of the first and second fragments.

# We are told one fragment fell near the gun. We can assume its landing distance is 0.
R1 = 0

# We can now set up the equation to find R2, the landing distance of the second fragment.
# 2 * I = (R1 + R2) / 2
# Rearranging for R2 gives: R2 = 4 * I - R1
# Since R1 = 0, the equation simplifies to: R2 = 4 * I

# Calculate the landing distance of the second fragment.
R2 = 4 * I

# Print the final calculation and the result.
# The question asks for the maximum distance one can be to remain safe,
# which corresponds to the landing spot of the second fragment.
print("The landing distance of the second fragment (R2) is found using the center of mass principle.")
print("The formula is: R2 = 4 * I")
print("Substituting the given value for I:")
print(f"R2 = 4 * {I} = {R2}")
print(f"\nTherefore, the second fragment will land at {R2} m from the gun. You must be beyond this distance to be safe.")
