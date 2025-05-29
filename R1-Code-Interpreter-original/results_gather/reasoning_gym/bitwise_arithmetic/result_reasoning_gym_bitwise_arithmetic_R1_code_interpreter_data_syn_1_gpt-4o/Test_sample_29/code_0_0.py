# Define the hexadecimal numbers
num1 = 0x4c17
num2 = 0x532a

# Perform the multiplication
product = num1 * num2

# Right shift the result by 2 bits (equivalent to integer division by 4)
result = product >> 2

# Print the result in hexadecimal format
print(hex(result))