# The problem is to find the maximum value of |b| + |c| for a quadratic
# P(x) = ax^2 + bx + c, given that |P(x)| <= 1 for all x in [-1, 1].

# Based on mathematical analysis, the maximum value is 1.5.
# This can be achieved with the polynomial P(x) = -0.5*x**2 + 1*x + 0.5.
# For this polynomial, the coefficients are:
a = -0.5
b = 1.0
c = 0.5

# The maximum value is |b| + |c|.
max_value = abs(b) + abs(c)

# We print the final equation as requested.
print("The maximum value of |b| + |c| is found to be 1.5.")
print("This is achieved, for instance, by the polynomial P(x) = -0.5x^2 + x + 0.5.")
print("The final equation demonstrating this maximum value is:")
print(f"|{b}| + |{c}| = {max_value}")
