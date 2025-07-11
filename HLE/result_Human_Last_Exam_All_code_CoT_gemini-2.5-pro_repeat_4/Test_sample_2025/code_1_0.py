import math

# The problem asks for the value of the probability density function f_Z(z) at z = 0.2.
# Based on the reasoning provided, the pdf of Z is given by the formula:
# f_Z(z) = -2 * (z*ln(z) + (1-z)*ln(1-z))

# We need to calculate this for z = 0.2.
z = 0.2

# Calculate the components of the formula.
val_z = z
val_1_minus_z = 1 - z
log_z = math.log(val_z)
log_1_minus_z = math.log(val_1_minus_z)

# Calculate the terms inside the parenthesis
term1 = val_z * log_z
term2 = val_1_minus_z * log_1_minus_z
total_term = term1 + term2

# Final result
result = -2 * total_term

# As requested, we show the numbers in the final equation.
print(f"The formula for the pdf is f(z) = -2 * (z*ln(z) + (1-z)*ln(1-z))")
print(f"We want to evaluate f(z) at z = {val_z}")
print(f"The equation becomes: f({val_z}) = -2 * ({val_z}*ln({val_z}) + {val_1_minus_z}*ln({val_1_minus_z}))")
print(f"Calculating the logarithms: ln({val_z}) = {log_z}, ln({val_1_minus_z}) = {log_1_minus_z}")
print(f"The equation with log values: f({val_z}) = -2 * ({val_z}*({log_z}) + {val_1_minus_z}*({log_1_minus_z}))")
print(f"Calculating the terms: {val_z}*ln({val_z}) = {term1}, {val_1_minus_z}*ln({val_1_minus_z}) = {term2}")
print(f"The equation becomes: f({val_z}) = -2 * ({term1} + {term2})")
print(f"Summing the terms inside the parenthesis: f({val_z}) = -2 * ({total_term})")
print(f"The final value is: f({val_z}) = {result}")
