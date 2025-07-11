import math

# We derived the maximal overhang to be (11 * sqrt(2)) / 12.
# The problem requires this to be in the form (a + sqrt(b)) / (1 + c).

# Let's convert the expression:
# Overhang = (11 * sqrt(2)) / 12
# We can write 11 * sqrt(2) as sqrt(11^2 * 2) = sqrt(121 * 2) = sqrt(242).
# So, Overhang = sqrt(242) / 12.

# Now we match this to the format (a + sqrt(b)) / (1 + c):
# Numerator: a + sqrt(b) = 0 + sqrt(242)  => a = 0, b = 242
# Denominator: 1 + c = 12                  => c = 11

a = 0
b = 242
c = 11

# The problem asks to output the integers a, b, c.
# The final equation is (0 + sqrt(242)) / (1 + 11).
print(f"{a} {b} {c}")