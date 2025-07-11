import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

z = 0.2

# The probability density function f_Z(z) for the random variable Z has been
# derived as f_Z(z) = 6*z - 6*z^2.
# We need to calculate the value of this function for z = 0.2.

# First term calculation
term1_coeff = 6
term1_val = term1_coeff * z

# Second term calculation
term2_coeff = 6
term2_val = term2_coeff * (z**2)

# Final result
result = term1_val - term2_val

print(f"The probability density function is f_Z(z) = 6*z - 6*z^2.")
print(f"We want to evaluate this at z = {z}.")
print(f"f_Z({z}) = {term1_coeff} * {z} - {term2_coeff} * {z}^2")
print(f"f_Z({z}) = {term1_val} - {term2_coeff} * {z*z}")
print(f"f_Z({z}) = {term1_val} - {term2_val}")
print(f"f_Z({z}) = {result}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)