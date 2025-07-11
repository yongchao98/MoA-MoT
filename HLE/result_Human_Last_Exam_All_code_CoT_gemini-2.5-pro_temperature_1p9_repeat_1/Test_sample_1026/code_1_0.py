import sys

# Define the given physical quantities
m = 0.20  # mass of the body in kg
M = 0.80  # mass of the guide in kg
R_cm = 20   # radius of the circular arcs in cm
d_cm = 50   # length of the straight section in cm

# Convert all length units to SI units (meters) for consistency
R = R_cm / 100.0
d = d_cm / 100.0

# The principle of conservation of the center of mass states that since there are no external
# horizontal forces, the displacement of the guide (ΔX_M) can be found using the formula:
# ΔX_M = - (m * Δx_rel) / (M + m)
# where Δx_rel is the horizontal displacement of the mass 'm' relative to the guide.

# The path consists of a quarter arc (horizontal span R), a straight section (d),
# and a second quarter arc (horizontal span R).
# So, the total horizontal displacement relative to the guide is Δx_rel = 2*R + d.
delta_x_rel = 2 * R + d

# Now, we calculate the numerator and denominator of the displacement formula.
numerator = -m * delta_x_rel
denominator = M + m

# Calculate the final displacement of the guide
delta_X_M = numerator / denominator

# Print the calculation steps with the numerical values
print("The formula for the guide's horizontal displacement (ΔX_M) is:")
print("ΔX_M = - (m * (2 * R + d)) / (M + m)\n")

print("Plugging in the values:")
print(f"ΔX_M = - ({m} kg * (2 * {R} m + {d} m)) / ({M} kg + {m} kg)")
print(f"ΔX_M = - ({m} * ({delta_x_rel})) / ({denominator})")
print(f"ΔX_M = {numerator} / {denominator}")
print(f"ΔX_M = {delta_X_M:.2f} m\n")

print(f"The horizontal displacement of the guide is {delta_X_M:.2f} meters (or {delta_X_M*100:.0f} cm to the left).")

# Output the final answer in the required format
# Redirecting final answer to stdout for the platform
sys.stdout.write(f"<<<{delta_X_M}>>>")