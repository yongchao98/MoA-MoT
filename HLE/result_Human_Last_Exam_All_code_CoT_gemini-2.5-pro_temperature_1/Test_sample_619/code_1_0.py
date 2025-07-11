import math

# The problem is to find a formula for the minimum value of the sum of sizes of the sets.
# Based on analysis, the formula that fits best for n>=2 is floor(n^2 / 4) + 2.
# This corresponds to answer choice B.
# The code below will print this formula. The problem does not specify a value for n.

# The final formula is floor(n^2 / 4) + 2
# The numbers in the equation are 2 and 4.
n_squared_str = "n^2"
denominator = 4
added_value = 2

# Printing the formula as requested.
print(f"The minimum value is given by the formula: floor({n_squared_str} / {denominator}) + {added_value}")
