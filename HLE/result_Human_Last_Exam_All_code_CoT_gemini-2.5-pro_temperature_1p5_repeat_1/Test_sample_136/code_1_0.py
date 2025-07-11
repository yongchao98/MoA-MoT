import math

# The derivation of the limit leads to the expression exp(-pi/2).
# Here is the Python code to compute this value.

# Define the constants in the equation
# pi is the mathematical constant, approximately 3.14159
# The expression is e raised to the power of (-pi / 2)
pi_val = math.pi
exponent_val = -pi_val / 2
final_result = math.exp(exponent_val)

# Output the steps of the final calculation
print("The limit is given by the expression: exp(-pi / 2)")
print(f"The value of pi is: {pi_val}")
print(f"The value of the exponent (-pi / 2) is: {exponent_val}")
print(f"The final result is exp({exponent_val}), which evaluates to:")
print(final_result)