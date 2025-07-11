# Define the given variables.
# I is the horizontal distance from the gun to the highest point of elevation in meters.
I = 500

# The problem states one fragment fell near the gun. To find the maximum
# distance for the second fragment, we assume the first fragment landed at the
# origin (x1 = 0).
x1 = 0

# The center of mass (CoM) of the system continues on the original trajectory.
# The total range of the original trajectory would have been 2 * I.
# This is the landing spot of the CoM.
R_CoM = 2 * I

# For two equal masses, the CoM is the average of their positions: R_CoM = (x1 + x2) / 2.
# We can rearrange this to solve for x2, the landing position of the second fragment:
# x2 = 2 * R_CoM - x1
# Substituting R_CoM = 2 * I and x1 = 0, we get:
# x2 = 2 * (2 * I) - 0 = 4 * I
max_distance = 4 * I

# Print the final equation with the numbers and the result.
print("The problem can be solved using the principle of conservation of the center of mass.")
print("The original trajectory's full range (and the landing spot of the Center of Mass) is R_CoM = 2 * I.")
print(f"R_CoM = 2 * {I} = {R_CoM} m.")
print("\nOne fragment lands at x1 = 0 m (near the gun).")
print("The landing spot of the Center of Mass is the average of the two fragments' landing spots: R_CoM = (x1 + x2) / 2.")
print("Solving for the second fragment's distance (x2): x2 = 2 * R_CoM - x1.")
print(f"x2 = 2 * {R_CoM} - {x1}")
print(f"x2 = {2 * R_CoM - x1} m.")
print("\nAlternatively, a simplified final equation is x2 = 4 * I.")
print(f"The maximum safe distance is calculated as: 4 * {I} = {max_distance} meters.")
