import numpy as np

# Define the polynomial for the intercepts
# P(x) = 3x^6 - 20x^5 - 2400x^4 + 1440000x^2 - 192000000
coeffs = [3, -20, -2400, 0, 1440000, 0, -192000000]

# Find the roots of the polynomial
roots = np.roots(coeffs)

# The intercepts are the real roots
real_roots = roots[np.isreal(roots)].real
x_max = max(abs(real_roots))
y_max = x_max # The polynomial for y-intercepts is the same

# The poles are at z = a + 2*pi*i*k
# a in {-2024, ..., 2024}
# k is an integer

# Count number of 'a' values inside the contour
# We need |a| < x_max
num_a = 2 * int(x_max) + 1

# Count number of 'k' values inside the contour
# We need |2*pi*k| < y_max  => |k| < y_max / (2*pi)
num_k = 2 * int(y_max / (2 * np.pi)) + 1

# Total number of poles inside the contour
total_poles = num_a * num_k

# The residue at each pole is 1.
# The integral is 2*pi*i * (sum of residues) = 2*pi*i * total_poles
integral_coeff = 2 * total_poles

print(f"The maximum x-intercept is at approximately: {x_max:.4f}")
print(f"The maximum y-intercept is at approximately: {y_max:.4f}")
print(f"Integer values of 'a' are in the range: -{int(x_max)} to {int(x_max)}")
print(f"Number of 'a' values: {num_a}")
print(f"The condition for k is |k| < {y_max / (2 * np.pi):.4f}")
print(f"Integer values of 'k' are in the range: -{int(y_max / (2 * np.pi))} to {int(y_max / (2 * np.pi))}")
print(f"Number of 'k' values: {num_k}")
print(f"Total number of poles enclosed by the contour Gamma is: {num_a} * {num_k} = {total_poles}")
print(f"The value of the contour integral is 2 * pi * i * (number of poles).")
print(f"The equation is: 2 * {total_poles} * pi * i")
print(f"So the result is {integral_coeff} * pi * i")
