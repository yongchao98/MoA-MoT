import math

# The maximal overhang is derived as (11 * sqrt(2)) / 12 side lengths.
# We need to express this in the form (a + sqrt(b)) / (1 + c).
# Overhang = (11 * sqrt(2)) / 12 = sqrt(121 * 2) / 12 = sqrt(242) / 12
# We can write this as (0 + sqrt(242)) / (1 + 11).

# From this form, we can identify the integer values for a, b, and c.
a = 0
b = 242
c = 11

# The problem asks to output the numbers for the final equation.
# The format from the example is "a b c".
print(f"{a} {b} {c}")

# This section is for verification and explanation, not part of the primary output.
# It demonstrates that the chosen values correctly represent the overhang.
# overhang_value = (11 * math.sqrt(2)) / 12
# formula_value = (a + math.sqrt(b)) / (1 + c)
# print(f"Calculated overhang: {overhang_value}")
# print(f"Value from formula: {formula_value}")
# print(f"Does sqrt(b) meet the condition? sqrt({b}) = {math.sqrt(b)}, which is non-integer: {math.sqrt(b) != int(math.sqrt(b))}")