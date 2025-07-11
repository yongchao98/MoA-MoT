import numpy as np

# Define the polynomial Q(x) with roots at 0, +/-1, +/-2, +/-3, +/-4
# Q(x) = x(x^2-1)(x^2-4)(x^2-9)(x^2-16)
# Expanding this gives:
# Q(x) = x^9 - 30x^7 + 273x^5 - 820x^3 + 576x
Q_coeffs = [1, 0, -30, 0, 273, 0, -820, 0, 576, 0]
Q = np.poly1d(Q_coeffs)

# Calculate the derivative Q'(x)
# Q'(x) = 9x^8 - 210x^6 + 1365x^4 - 2460x^2 + 576
Q_prime = np.polyder(Q)

# To find the extrema of Q'(x), we find the roots of its derivative Q''(x)
Q_double_prime = np.polyder(Q_prime)
critical_points = np.roots(Q_double_prime)

# Evaluate Q'(x) at these critical points and endpoints if considering an interval.
# Since Q' is a polynomial of even degree with positive leading coeff, it has a global minimum.
# We are interested in the global minimum and maximum of Q'(x).
# Real critical points for Q'(x) are the real roots of Q''(x).
real_critical_points = critical_points[np.isreal(critical_points)].real
extrema_values = Q_prime(real_critical_points)

# As x -> +/- infinity, Q'(x) -> +infinity, so there is no global maximum.
# Let's find the minimum value.
min_Q_prime = np.min(extrema_values)

# The polynomial Q(x) is odd, so Q'(x) is even.
# We can analyze Q'(x) by substituting u = x^2.
# Q'(x) = 9u^4 - 210u^3 + 1365u^2 - 2460u + 576 for u >= 0
# A plot would show Q'(x) has local maxima and minima.
# For our purpose, let's find the absolute min and max of Q'(x) on the interval containing the roots, e.g., [-4, 4].
# A more robust numerical approach:
x_vals = np.linspace(-4, 4, 10000)
q_prime_vals = Q_prime(x_vals)
min_val = np.min(q_prime_vals)
max_val = np.max(q_prime_vals)

# We need c * Q'(x) > -1.
# This requires two conditions if c is positive:
# 1) c * min_val > -1  => c < -1/min_val
# 2) We also need to handle the maximum value if we want to bound h'(x).
# However, for c > 0, c * max_val will be positive, so we only need to worry about the minimum.
# If c is negative:
# 1) c * min_val > -1 => c > -1/min_val (since min_val is negative, -1/min_val is positive)
# 2) c * max_val > -1 => c > -1/max_val (since max_val is positive, -1/max_val is negative)
# So we need c to be in (-1/max_val, 0)

c = -1 / (max_val * 2) # A safe choice for a negative c

# Let's verify the condition h'(x) = 1 + c*Q'(x) > 0
# min(h'(x)) = 1 + c*max(Q'(x)) = 1 + c*max_val = 1 + (-1/(2*max_val))*max_val = 1 - 0.5 = 0.5 > 0
# max(h'(x)) = 1 + c*min(Q'(x)) = 1 + (-1/(2*max_val))*min_val. Since min_val is negative, this is > 1.

h_coeffs = np.polyadd([0, 0, 0, 0, 0, 0, 0, 0, 1, 0], c * Q_coeffs)

print("We need to find the maximum number of solutions to f(g(x)) = x.")
print("This is equivalent to finding the roots of h(x) - x = 0, where h(x) is a degree 9 polynomial.")
print("We can construct such a polynomial h(x) = x + c*Q(x) where Q(x) has 9 real roots.")
print("Let Q(x) have roots at 0, +/-1, +/-2, +/-3, +/-4.")
print(f"Q(x) = {Q}")
print(f"Its derivative is Q'(x) = {Q_prime}")
print(f"The minimum value of Q'(x) on [-4, 4] is approximately {min_val:.2f}.")
print(f"The maximum value of Q'(x) on [-4, 4] is approximately {max_val:.2f}.")
print(f"We can choose a negative c, for example c = {c:.6f}.")
print("With this c, h(x) = x + c*Q(x) has 9 fixed points, and its derivative h'(x) = 1 + c*Q'(x) is always positive.")
print("The final constructed polynomial with 9 fixed points is h(x) - x = c*Q(x) = 0.")
final_poly_coeffs = c * Q_coeffs
final_poly = np.poly1d(final_poly_coeffs)
# We need to print the equation. f(g(x)) = x is h(x)=x. So h(x)-x=0.
print("This leads to the equation having 9 roots:")
print(f"{final_poly} = 0")
print("\nThus, it is possible to have 9 fixed points.")
