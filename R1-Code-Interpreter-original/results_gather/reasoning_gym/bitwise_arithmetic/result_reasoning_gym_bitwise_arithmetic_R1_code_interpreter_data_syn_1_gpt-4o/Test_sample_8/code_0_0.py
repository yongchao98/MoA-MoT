# Define the hexadecimal numbers
num1 = 0xab20
num2 = 0xf737

# Perform the multiplication
product = num1 * num2

# Right shift the result by 2 bits (equivalent to integer division by 4)
result = product >> 2

# Print the result in hexadecimal format
# If the result is negative, it will be automatically represented as a negative value
print(hex(result))