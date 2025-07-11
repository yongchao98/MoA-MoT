import math

# Define constants
pi = math.pi
gamma = 0.57721566490153286060651209008240243104215933593992

# Parameters from the problem
dist_z = 3600.0

# Calculate the terms in the denominator of the probability formula
log_dist = math.log(dist_z)
half_log_8 = math.log(8) / 2
denominator = 4 * (log_dist + gamma + half_log_8)

# Calculate the final probability
probability = pi / denominator

print("The probability is calculated using the formula:")
print("p = pi / (4 * (log(||z||) + gamma + log(8)/2))")
print("\nWith the following values:")
print(f"pi = {pi:.4f}")
print(f"||z|| = {dist_z:.0f}")
print(f"log(||z||) = {log_dist:.4f}")
print(f"gamma = {gamma:.4f}")
print(f"log(8)/2 = {half_log_8:.4f}")

print("\nThe final equation with numbers is:")
print(f"p = {pi:.4f} / (4 * ({log_dist:.4f} + {gamma:.4f} + {half_log_8:.4f}))")
print(f"p = {pi:.4f} / (4 * {log_dist + gamma + half_log_8:.4f})")
print(f"p = {pi:.4f} / {denominator:.4f}")
print(f"\nThe calculated probability is: {probability:.4f}")
print(f"The approximate answer with two significant digits is: {probability:.2g}")
print(f"<<<{probability:.2g}>>>")