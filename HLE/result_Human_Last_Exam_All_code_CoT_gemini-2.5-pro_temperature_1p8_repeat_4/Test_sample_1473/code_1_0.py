import math

# The value of the integral is given by the expression I = pi * ln(1 + sqrt(2)).
# Here, we calculate this value.

# Define the constants from the equation I = pi * ln(1 + sqrt(2))
pi_val = math.pi
num_1 = 1
num_2 = 2

# Calculate the components of the expression
sqrt2_val = math.sqrt(num_2)
inner_term = num_1 + sqrt2_val
log_val = math.log(inner_term)

# Calculate the final value of the integral I
I = pi_val * log_val

# Print the components and the final result based on the equation I = pi * ln(1 + sqrt(2))
print("The final equation is: I = pi * ln(1 + sqrt(2))")
print(f"Value of pi: {pi_val}")
print(f"Value of 1: {num_1}")
print(f"Value of 2: {num_2}")
print(f"Value of sqrt(2): {sqrt2_val}")
print(f"Value of 1 + sqrt(2): {inner_term}")
print(f"Value of ln(1 + sqrt(2)): {log_val}")
print("---")
print(f"Final Value of I = {pi_val} * {log_val}")
print(f"I = {I}")
