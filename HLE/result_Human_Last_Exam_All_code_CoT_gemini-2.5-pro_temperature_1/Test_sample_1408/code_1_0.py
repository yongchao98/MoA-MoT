import math

# Based on the derivation, the maximal overhang is 11*sqrt(2)/12.
# We are asked to represent this in the format (a + sqrt(b)) / (1 + c).

# We rewrite the overhang as:
# (11 * sqrt(2)) / 12 = sqrt(121 * 2) / 12 = sqrt(242) / 12
#
# Comparing this with the target format (a + sqrt(b)) / (1 + c), we find:
# a = 0
# b = 242
# 1 + c = 12  =>  c = 11

a = 0
b = 242
c = 11

# The problem asks to output each number in the final equation.
print(f"The maximal overhang is ({a} + \u221A{b}) / (1 + {c}) which is equal to {11 * math.sqrt(2) / 12:.5f}")
# The result in the format "a b c"
print(f"{a} {b} {c}")
