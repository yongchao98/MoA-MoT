import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# Value of z at which to evaluate the PDF
z = 0.2

# The PDF is f(z) = 6 * z * (1 - z)
# Coefficients from the derived formula
a = 6
b = 1

# Calculate the PDF value
result = a * z * (b - z)

# Print the equation with the numbers plugged in, as requested
print(f"{a} * {z} * ({b} - {z}) = {result}")

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)

# Final answer format
final_answer = f"<<<{result}>>>"