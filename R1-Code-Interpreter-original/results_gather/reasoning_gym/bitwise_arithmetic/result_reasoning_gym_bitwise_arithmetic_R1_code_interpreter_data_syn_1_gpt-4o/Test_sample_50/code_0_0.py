# Define the hexadecimal numbers
num1 = 0xf906
num2 = 0xd416

# Multiply the numbers
product = num1 * num2

# Left shift the result by 1
result = product << 1

# Convert the result to hexadecimal
hex_result = hex(result)

# Check if the result is negative
# Assuming a 32-bit signed integer for this example
if result >= 2**31:
    # Convert to negative value
    result = result - 2**32
    hex_result = hex(result)

print(hex_result)