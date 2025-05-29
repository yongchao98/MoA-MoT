# Step 1: Convert base-15 number '2c1' to decimal
base_15_number = '2c1'
base_15_digits = '0123456789abcde'
decimal_value = 0

# Calculate the decimal value
for i, digit in enumerate(reversed(base_15_number)):
    decimal_value += base_15_digits.index(digit) * (15 ** i)

# Step 2: Convert the decimal value to hexadecimal
hexadecimal_value = ''
hex_digits = '0123456789abcdef'

while decimal_value > 0:
    remainder = decimal_value % 16
    hexadecimal_value = hex_digits[remainder] + hexadecimal_value
    decimal_value //= 16

print(hexadecimal_value)