# The maximal overhang for three identical cubes is 1 block sidelength.
# This is achieved by a non-linear stacking arrangement where one cube acts as a counterweight.
#
# Configuration:
# - Block C2: CM at x=0, on the table edge.
# - Block C3: CM at x=-1/2, on the table, as a counterweight.
# - Block C1: CM at x=1/2, stacked on top of C2.
#
# Stability analysis:
# - CM of C1 is on the edge of C2's support: Stable.
# - Combined CM of {C1, C2} acts at x=1/4, which is within C2's base on the table: Stable.
# - Combined CM of {C1, C2, C3} is at x=0, which is on the table edge: Stable.
#
# The rightmost point is the right edge of C1, at x = 1/2 (CM) + 1/2 (half length) = 1.
# So, the maximal overhang is 1.
#
# We need to express 1 in the format (a + sqrt(b)) / (1 + c).
# 1 = (1 + sqrt(0)) / (1 + 0)
#
# From this, we identify the integers:
# a = 1
# b = 0
# c = 0
#
# Checking the constraints:
# - a, b, c are non-negative integers (1, 0, 0 are).
# - sqrt(b) is zero or noninteger (sqrt(0) = 0 is zero).
# - c is minimal (c=0 is the minimal non-negative integer).

a = 1
b = 0
c = 0

print(f"The equation for the maximal overhang is ({a} + sqrt({b})) / (1 + {c})")
print(a, b, c)
