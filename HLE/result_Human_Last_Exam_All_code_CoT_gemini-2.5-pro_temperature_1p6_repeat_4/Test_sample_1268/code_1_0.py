import math

# The problem asks for the upper bound for a norm in a number field k,
# constructed from a squarefree natural number N, in relation to the covolume V.
# This is interpreted as the Minkowski bound for a quadratic field.

# The Minkowski bound, B, is given by the inequality B <= C * V.
# The constant C depends on whether the quadratic field is real or imaginary.
# - For real quadratic fields, C = 1/2.
# - For imaginary quadratic fields, C = 4/pi.

# To find a single, universal upper bound that applies to all quadratic fields,
# we must take the larger of the two constants.
constant_C = max(0.5, 4 / math.pi)

# The numbers that form this constant are 4 and pi.
numerator = 4
denominator_value = math.pi

# The final code outputs the equation for the upper bound,
# explicitly showing the numbers involved in the calculation.
print("The universal upper bound (B) for the ideal norm in any quadratic field is related to the covolume (V) by the following inequality:")

# The final equation is: B <= (4 / pi) * V
print("\nFinal Equation:")
print(f"B <= ({numerator} / {denominator_value}) * V")

# We also show the calculated value of the constant.
print(f"\nWhich simplifies to:")
print(f"B <= {constant_C} * V")