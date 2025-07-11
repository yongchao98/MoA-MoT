import numpy as np

# Step 1: Determine the azimuthal winding number P.
# The field's azimuthal angle is f = atan2(y, x), which is the same as the
# cylindrical coordinate phi. As we circle the z-axis, f sweeps 2*pi once.
P = 1

# Step 2: Determine the polar winding number W.
# This depends on the range of the polar angle G.
# G = PI * exp(-10 * r2)
# We need to find the minimum and maximum of G.

# G is maximum when r2 is minimum.
# r2 = sqrt((x*x+y*y-0.5)*(x*x+y*y-0.5)+z*z)
# The minimum of r2 is 0 (at the circle x^2+y^2=0.5, z=0).
r2_min = 0
G_max = np.pi * np.exp(-10 * r2_min)

# G is minimum when r2 is maximum.
# The maximum of r2 is infinity (at spatial infinity).
# r2_max -> infinity
G_min = 0.0

# Calculate W using the formula W = |G_max - G_min| / PI.
W = abs(G_max - G_min) / np.pi

# Step 3: Calculate the total Hopf charge Q = P * W.
Q = P * W

# Print the results following the formula Q = P * W
print("The Hopf charge Q is calculated as Q = P * W.")
print(f"Azimuthal winding number, P = {P}")
print(f"Polar angle G at the structure's core (r2=0), G_max = {G_max:.4f}")
print(f"Polar angle G at spatial infinity (r2=inf), G_min = {G_min:.4f}")
print(f"Polar winding number, W = (|G_max - G_min|) / PI = {W}")
print("Final Equation:")
print(f"{Q:.0f} = {P} * {W:.0f}")
