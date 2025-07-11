import math

# Step 1: Define the known values from Frame 2 based on the problem statement.
# Angle between S1 and S2 is pi/2.
cos_theta_prime_12 = math.cos(math.pi / 2) 
# Angle between S2 and S3 is 3*pi/4.
cos_theta_prime_23 = math.cos(3 * math.pi / 4)

# Step 2: Use the Lorentz invariant cross-ratio.
# In Frame 1, the cross-ratio is 1. So it must be 1 in Frame 2.
# (1-cos(theta'_14))(1-cos(theta'_23)) = (1-cos(theta'_12))(1-cos(theta'_34))
# We want to find the ratio R = (1-cos(theta'_14)) / (1-cos(theta'_34)).
# Rearranging the equation gives: R = (1-cos(theta'_12)) / (1-cos(theta'_23))

# Step 3: Calculate the numerator and denominator of the ratio.
numerator_val = 1 - cos_theta_prime_12
denominator_val = 1 - cos_theta_prime_23
result = numerator_val / denominator_val

# Step 4: Print the final equation with the computed values.
# Note: cos(pi/2) is exactly 0.
print("The problem reduces to calculating the ratio R = (1 - cos(theta'_12)) / (1 - cos(theta'_23))")
print("Using the given angles in the second frame:")
print(f"cos(theta'_12) = cos(pi/2) = {cos_theta_prime_12:.4f}")
print(f"cos(theta'_23) = cos(3*pi/4) = {cos_theta_prime_23:.4f}")
print("\nThe final equation with plugged-in values is:")
print(f"R = (1 - {cos_theta_prime_12:.1f}) / (1 - ({cos_theta_prime_23:.4f}))")
print(f"R = {numerator_val:.1f} / {denominator_val:.4f}")
print(f"R = {result:.4f}")
print("\nThe exact value is 2 - sqrt(2).")
