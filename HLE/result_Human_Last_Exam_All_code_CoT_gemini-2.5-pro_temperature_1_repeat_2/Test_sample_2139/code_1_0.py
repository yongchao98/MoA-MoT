import math

# Given initial value y(0)
# y(0) = (128 * 3**(1/6) * Gamma(2/3))**-1
# We can compute this constant first.
gamma_2_3 = math.gamma(2/3)
y0_inv = 128 * (3**(1/6)) * gamma_2_3
y0 = 1 / y0_inv

# The hypothesized solution form is y(t) = y(0) / (t**2 * sec(t) + cos(t))
# We need to calculate y(pi/4)

t = math.pi / 4

# Denominator at t = pi/4
# f(pi/4) = (pi/4)**2 * sec(pi/4) + cos(pi/4)
# sec(pi/4) = sqrt(2), cos(pi/4) = 1/sqrt(2)
sec_pi_4 = math.sqrt(2)
cos_pi_4 = 1 / math.sqrt(2)
denominator = (t**2) * sec_pi_4 + cos_pi_4

# y(pi/4) = y(0) / denominator
y_pi_4 = y0 / denominator

# The full expression is y(pi/4) = (8 * sqrt(2) * y(0)) / (pi**2 + 8)
# Let's calculate using this simplified symbolic form to be sure.
pi_sq = math.pi**2
y_pi_4_symbolic = (8 * math.sqrt(2) * y0) / (pi_sq + 8)

# The results from both methods should be identical.
# Let's print the result and the equation leading to it.

print(f"The initial value y(0) is approximately: {y0:.6f}")
print(f"The value of the denominator f(pi/4) is: (pi/4)^2 * sec(pi/4) + cos(pi/4) = {denominator:.6f}")
print(f"The radius of the balloon at t=pi/4 is y(pi/4) = y(0) / f(pi/4)")
print(f"y(pi/4) = {y0:.6f} / {denominator:.6f}")
print(f"The final calculated value is:")
print(f"{y_pi_4}")
