import math

# Define the constants and given values for the equation
# The simplified ratio of the differential cross-sections is Ratio = 3 / sin^2(theta)
# The constant numerator
c_num = 3
# The scattering angle in radians
theta = math.pi / 30

# Calculate sin(theta)
sin_theta = math.sin(theta)

# Calculate the final ratio
ratio = c_num / (sin_theta ** 2)

# Print the equation with the plugged-in numbers
print("The ratio of the differential cross-sections is given by the formula:")
print("Ratio = (Monopole Cross-Section) / (Dipole Cross-Section)")
print("After substitution and simplification, this becomes:")
print(f"Ratio = 3 / sin^2({theta:.5f})")
print(f"Ratio = 3 / ({sin_theta:.5f})^2")
print(f"Ratio = 3 / {sin_theta**2:.5f}")
print(f"Final calculated ratio is:")
print(ratio)