import math

# Given inner products (cosines of the angles)
c_hb = 0.9375  # lim <h_p, b_p>
c_hz = 0.9     # lim <h_p, z_p>

# The leading eigenvector h_p is expected to lie between the signal vectors b_p and z_p.
# This implies the angle between b_p and z_p is the sum of the other two angles:
# theta_bz = theta_hb + theta_hz
#
# Using the cosine of a sum formula:
# cos(theta_bz) = cos(theta_hb) * cos(theta_hz) - sin(theta_hb) * sin(theta_hz)
#
# We can find the sines from the cosines: sin(theta) = sqrt(1 - cos^2(theta))
# (sine is non-negative as angles are in [0, pi])

s_hb = math.sqrt(1 - c_hb**2)
s_hz = math.sqrt(1 - c_hz**2)

# Calculate the final cosine value
c_bz = c_hb * c_hz - s_hb * s_hz

print("The equation we are solving is:")
print(f"cos(theta_bz) = cos(theta_hb) * cos(theta_hz) - sin(theta_hb) * sin(theta_hz)")
print(f"cos(theta_bz) = {c_hb} * {c_hz} - sqrt(1 - {c_hb}^2) * sqrt(1 - {c_hz}^2)")
print(f"cos(theta_bz) = {c_hb * c_hz} - {s_hb} * {s_hz}")
print("\nThe final result is:")
print(c_bz)