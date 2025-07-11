import math

# Based on the analysis, the volume required is half the volume of the ellipsoid S2,
# because the solid ellipsoid is fully contained within the tangent cone.
# The volume of interest is for y > 0.

# Define the parameters of the ellipsoid S2: x^2/3 + y^2/4 + z^2/3 = 1
a_sq = 3
b_val = 2
c_sq = 3

# String representations for printing
a_str = f"sqrt({a_sq})"
c_str = f"sqrt({c_sq})"

# Calculate the final volume.
# V = (1/2) * V_ellipsoid = (1/2) * (4/3) * pi * a * b * c
# V = (1/2) * (4/3) * pi * sqrt(3) * 2 * sqrt(3) = (1/2) * (4/3) * pi * 6 = 4 * pi
final_coefficient = 4
final_volume = final_coefficient * math.pi

# Print the explanation and the final calculation as requested.
print("The volume of the space is the volume of the part of the ellipsoid with y > 0.")
print("The semi-axes of the ellipsoid are a = sqrt(3), b = 2, and c = sqrt(3).")
print("The volume of the full ellipsoid is V_ellipsoid = (4/3) * pi * a * b * c.")
print("The requested volume V is half of the ellipsoid's total volume.")

print("\nThe final equation with all numbers is:")
print(f"V = (1/2) * (4/3) * \u03C0 * {a_str} * {b_val} * {c_str}")

# Show the simplification step
product_abc = 6
print(f"V = (1/2) * (4/3) * \u03C0 * {product_abc}")

# Show the final result in terms of pi
print(f"V = {final_coefficient} * \u03C0")

# Show the final numerical result
print(f"The numerical value is approximately {final_volume:.4f}")