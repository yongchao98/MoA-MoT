import math

# --- Setup ---
# The volume of the playdough is given.
volume_V = 1.0  # in cubic meters

# --- Derivation ---
# To maximize gravity at a point A (the origin), the optimal shape is bounded by the surface r = k * sqrt(cos(theta)).
# The volume of this shape is derived from integration: V = (4 * pi * k^3) / 15.
# We are given V = 1, so we can solve for k.
# 1 = (4 * pi * k^3) / 15
# k^3 = 15 / (4 * pi)
# The furthest distance from A (the origin) to the surface is the maximum value of r.
# r is maximized when cos(theta) is maximum (which is 1). So, r_max = k.
# Therefore, the furthest distance is k = (15 / (4 * pi))^(1/3).

# --- Calculation ---
pi_val = math.pi
k_cubed = (15 * volume_V) / (4 * pi_val)
r_max = k_cubed**(1/3)

# --- Output ---
# As requested, output the final equation with numerical values.
print("To find the furthest distance, we first solve for the constant 'k' which defines the optimal shape.")
print(f"The equation for the volume (V) of the shape is: V = (4 * pi * k^3) / 15")
print(f"Given V = {volume_V:.1f} m^3, we solve for k:")
print(f"1.0 = (4 * {pi_val:.6f} * k^3) / 15")
print(f"k^3 = (15 * 1.0) / (4 * {pi_val:.6f})")
print(f"k^3 = 15 / {4 * pi_val:.6f} = {k_cubed:.6f}")
print("")
print("The furthest distance (r_max) from point A to the surface is equal to k:")
print(f"r_max = k = ({k_cubed:.6f})^(1/3)")
print("\nThe final answer in meters is:")
print(r_max)
