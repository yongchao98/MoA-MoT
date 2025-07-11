import math

# Define the given horizontal distance to the highest point of elevation.
I = 500  # in meters

# --- Step 1: Determine the landing point of the Center of Mass (CM) ---
# The Center of Mass continues on the original trajectory.
# The total range of the original trajectory would be twice the distance to the apex.
R_cm = 2 * I

# --- Step 2: Define the landing point of the first fragment ---
# One fragment "fell near the gun", meaning it landed at the origin.
R1 = 0  # in meters

# --- Step 3: Use the Center of Mass formula to find the second fragment's position ---
# For two fragments of equal mass, the CM position is the average of their positions:
# R_cm = (R1 + R2) / 2
# We can rearrange this to solve for R2, the landing position of the second fragment:
# 2 * R_cm = R1 + R2
# R2 = 2 * R_cm - R1

R2 = 2 * R_cm - R1

# --- Step 4: Print the final result and the equation ---
# The maximum distance from the gun to be safe is the landing spot of the second fragment.
print("This problem can be solved using the principle of the center of mass.")
print(f"The horizontal distance to the apex is I = {I} m.")
print(f"The center of mass of the system will land at a distance of 2 * I = {2 * I} m from the gun.")
print(f"The first fragment lands at R1 = {R1} m from the gun.")
print("The landing position of the second fragment (R2) can be found using the equation: R2 = 2 * R_cm - R1.")
print(f"Substituting the values: R2 = 2 * {R_cm} - {R1}")
print(f"The final equation for the landing distance of the second fragment is: {R2} = 4 * {I}")

print("\n--- Final Answer ---")
print(f"The maximum distance from the gun where the second fragment can land is {R2} m.")
