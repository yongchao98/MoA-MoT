import numpy as np

# Coefficients of the quartic polynomial in u, from highest to lowest power.
# The equation is 49*u^4 - 392*u^3 + 1170*u^2 - 1544*u + 760 = 0
coeffs = [49, -392, 1170, -1544, 760]

print("The final polynomial equation to solve is:")
print(f"{coeffs[0]}*u^4 + ({coeffs[1]})*u^3 + {coeffs[2]}*u^2 + ({coeffs[3]})*u + {coeffs[4]} = 0")

# Find the roots of the polynomial
roots = np.roots(coeffs)

# Filter for real roots (roots where the imaginary part is close to zero)
# and count the number of unique real roots.
real_roots = []
tolerance = 1e-9
for root in roots:
    if abs(root.imag) < tolerance:
        real_part = root.real
        # Check if this root is already in our list (within tolerance)
        is_duplicate = False
        for r in real_roots:
            if abs(r - real_part) < tolerance:
                is_duplicate = True
                break
        if not is_duplicate:
            real_roots.append(real_part)

num_intersection_points = len(real_roots)

print(f"\nThe number of distinct real roots for u is: {num_intersection_points}")
print("This corresponds to the number of real intersection points between the two quadrics.")
print(f"Number of real intersection points: {num_intersection_points}")
