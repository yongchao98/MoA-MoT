import numpy as np

# Part 1: Find the root c = cos(theta) for the imaginary axis crossing
# The equation is 3*c^4 - 8*c^3 + 6*c^2 - 1 = 0
poly_coeffs = [3, -8, 6, 0, -1]
roots = np.roots(poly_coeffs)

# Filter for the real root c in [-1, 1] that is not 1.
c = None
for r in roots:
    # Check if the root is real and within the valid range for cos(theta)
    if np.isreal(r) and -1 <= r.real <= 1:
        # Exclude the trivial root at c=1 (theta=0)
        if not np.isclose(r.real, 1.0):
            c = r.real
            break

print(f"The equation for the imaginary axis crossing is 3*c^4 - 8*c^3 + 6*c^2 - 1 = 0, where c = cos(theta).")
print(f"A non-trivial real root for c is approximately: {c:.6f}")
print(f"The exact value of this root is -1/3.")

# Part 2: Calculate the angle based on this root
# We hypothesize the answer is related to tan(theta)
# tan(theta) = sin(theta) / cos(theta) = sqrt(1-c^2) / c
# Since cos(theta) is negative and theta for stability is in (pi/2, pi), tan(theta) is negative.
# We are interested in the positive angle alpha.

val_sqrt_8 = 2 * np.sqrt(2)
alpha_rad = np.arctan(val_sqrt_8)

print("\nFrom c = -1/3, we can calculate tan(theta) = -sqrt(8) = -2*sqrt(2).")
print("The stability angle alpha can be related to this value.")
print(f"We propose the exact value is arctan(2*sqrt(2)).")

# Final exact value expression
# Let's show the final output as requested.
# The angle is alpha, and we need to output the equation for it.
print("\nThe exact value of the angle alpha in terms of arctan() is given by the equation:")
print("alpha = arctan(2 * sqrt(2))")
print("\nNumerical value:")
print(f"2 * sqrt(2) is approximately {val_sqrt_8:.6f}")
print(f"alpha is approximately {alpha_rad:.6f} radians or {np.rad2deg(alpha_rad):.6f} degrees.")