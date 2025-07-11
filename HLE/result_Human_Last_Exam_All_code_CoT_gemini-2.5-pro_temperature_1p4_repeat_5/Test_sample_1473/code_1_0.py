import math

# The value of the integral I is given by the equation:
# I = pi * ln(1 + sqrt(2))

# Let's calculate the numerical value of the components and the final result.

# The numbers in the final equation
pi_val = math.pi
sqrt2_val = math.sqrt(2)
inner_term = 1 + sqrt2_val

# Final result
result = pi_val * math.log(inner_term)

print("The final equation for the integral is: I = pi * ln(1 + sqrt(2))")
print(f"Value of pi: {pi_val}")
print(f"Value of sqrt(2): {sqrt2_val}")
print(f"Value of 1 + sqrt(2): {inner_term}")
print("-" * 30)
print(f"The numerical value of the integral is: {result}")