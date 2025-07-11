import math

# Given parameters
m = 0.1  # kg
M = 10.0   # kg
theta_deg = 30.0
h = 2.0    # m
mu = 1.0 / 2.0
v0 = 1.0   # m/s (This information is not needed for the calculation of time)
g = 10.0   # m/s^2

# Convert angle to radians for math functions
theta = math.radians(theta_deg)
sin_theta = math.sin(theta)
cos_theta = math.cos(theta)

print("Step 1: Determine the relative acceleration of the block ('a_rel').")
print("The formula for a_rel, derived from Newton's laws and kinematic constraints, is:")
print("a_rel = g / (sin(theta) + (M / (M + m)) * cos(theta) * (cos(theta) + mu * sin(theta)) / (sin(theta) + mu * cos(theta)))")

# Calculate the terms for clarity
term1 = sin_theta
term_frac_num = cos_theta * (cos_theta + mu * sin_theta)
term_frac_den = sin_theta + mu * cos_theta
term2 = (M / (M + m)) * term_frac_num / term_frac_den
denominator = term1 + term2

# Calculate a_rel
a_rel = g / denominator

print("\n--- Calculation Details for a_rel ---")
print(f"theta = {theta_deg} degrees")
print(f"sin(theta) = {sin_theta:.4f}")
print(f"cos(theta) = {cos_theta:.4f}")
print(f"M / (M + m) = {M} / ({M} + {m}) = {M / (M + m):.4f}")
print("Numerator of the fraction term = cos(theta)*(cos(theta) + mu*sin(theta))")
print(f" = {cos_theta:.4f} * ({cos_theta:.4f} + {mu}*{sin_theta:.4f}) = {term_frac_num:.4f}")
print("Denominator of the fraction term = sin(theta) + mu*cos(theta)")
print(f" = {sin_theta:.4f} + {mu}*{cos_theta:.4f} = {term_frac_den:.4f}")
print("The complete denominator of the a_rel expression = sin(theta) + ...")
print(f" = {term1:.4f} + {M / (M + m):.4f} * {term_frac_num:.4f} / {term_frac_den:.4f} = {denominator:.4f}")

print(f"\na_rel = {g} / {denominator:.4f} = {a_rel:.4f} m/s^2")

print("\nStep 2: Calculate the distance 'd' the block slides down the incline.")
# d = h / sin(theta)
d = h / sin_theta
print(f"d = h / sin(theta) = {h} / {sin_theta:.4f} = {d:.4f} m")


print("\nStep 3: Calculate the time 't' using kinematics.")
print("The block starts from rest relative to the wedge, so we use d = (1/2) * a_rel * t^2.")
print("Solving for t: t = sqrt(2 * d / a_rel)")

# t = sqrt(2 * d / a_rel)
time_squared = 2 * d / a_rel
t = math.sqrt(time_squared)

print("\n--- Final Calculation for time 't' ---")
print(f"t = sqrt(2 * {d:.4f} / {a_rel:.4f})")
print(f"t = sqrt({2 * d:.4f} / {a_rel:.4f})")
print(f"t = sqrt({time_squared:.4f})")
print(f"t = {t:.4f} s")