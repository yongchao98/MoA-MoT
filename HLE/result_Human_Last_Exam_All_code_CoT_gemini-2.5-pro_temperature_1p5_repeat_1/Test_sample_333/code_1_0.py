# Define the given variables
# I is the horizontal distance from the gun to the highest point of elevation in meters.
I = 500
# x1 is the landing position of the first fragment in meters.
x1 = 0

print("Step 1: The center of mass (CM) of the system continues on the original trajectory.")
print(f"The horizontal distance to the highest point is given as I = {I} m.")
print("")

print("Step 2: Calculate the total horizontal range the CM would travel.")
# For a symmetric projectile path, the total range is twice the distance to the peak.
R_CM = 2 * I
print(f"The total range of the CM is R_CM = 2 * I = 2 * {I} = {R_CM} m.")
print("")

print("Step 3: Use the center of mass formula to find the landing spot of the second fragment.")
print("The projectile splits into two equal parts, so their masses are equal.")
print("The landing position of the first fragment is near the gun, so x1 = 0 m.")
print("The formula for the CM's landing spot is: R_CM = (x1 + x2) / 2")
print("We can solve for x2, the landing position of the second fragment: x2 = 2 * R_CM - x1")
print("")

print("Step 4: Substitute the values and calculate the final distance.")
# Calculate the landing position of the second fragment
x2 = 2 * R_CM - x1

# Print the final calculation and the result
print("The calculation is:")
print(f"x2 = 2 * {R_CM} - {x1}")
print(f"x2 = {x2}")
print("")
print(f"The maximum distance from the gun to be safe from the second fragment is {x2} meters.")
