import math

# Plan:
# 1. Define the constants needed for the calculation.
# 2. Calculate the distance of the starting point from the origin.
# 3. Apply the derived formula for the probability.
# 4. Print the steps of the calculation and the final result.

# Constants
gamma = math.gamma  # This is the Gamma function, we need the Euler-Mascheroni constant
euler_mascheroni_gamma = 0.5772156649
ln_8 = math.log(8)
start_pos_x = 3000
start_pos_y = 4000

# Calculate the distance |x_0|
r = math.sqrt(start_pos_x**2 + start_pos_y**2)
ln_r = math.log(r)

# Calculate the constant C = gamma + ln(8)
C = euler_mascheroni_gamma + ln_8

# Calculate the final probability using the formula P = ln(r) / (ln(r) + C)
probability = ln_r / (ln_r + C)

# Output the explanation and the result
print("The problem is to find the probability that a 2D random walk, conditioned to avoid the origin,")
print("starting at (3000, 4000), will never visit the four neighbors of the origin.")
print("\nThis probability P can be approximated by the formula:")
print("P = ln(r) / (ln(r) + C)")
print("\nWhere:")
print(f"r is the distance from the start position to the origin.")
print(f"C is a constant from the potential kernel theory, C = gamma + ln(8).")

print("\nFirst, we calculate the values of the components in the formula:")
print(f"The distance r = sqrt({start_pos_x}^2 + {start_pos_y}^2) = {r}")
print(f"The natural logarithm of r, ln({r:.0f}), is approximately: {ln_r:.4f}")
print(f"The Euler-Mascheroni constant, gamma, is approximately: {euler_mascheroni_gamma:.4f}")
print(f"The natural logarithm of 8, ln(8), is approximately: {ln_8:.4f}")
print(f"The constant C = {euler_mascheroni_gamma:.4f} + {ln_8:.4f} is approximately: {C:.4f}")

print("\nNow, we plug these numbers into the final equation:")
print(f"P = {ln_r:.4f} / ({ln_r:.4f} + {C:.4f})")
print(f"P = {ln_r:.4f} / {ln_r + C:.4f}")
print(f"P is approximately: {probability:.4f}")

print(f"\nThe approximate answer with two significant digits is: {probability:.2f}")
