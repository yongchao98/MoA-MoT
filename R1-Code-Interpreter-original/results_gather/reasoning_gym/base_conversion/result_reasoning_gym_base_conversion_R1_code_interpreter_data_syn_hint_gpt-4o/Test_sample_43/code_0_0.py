# Step 1: Convert base-11 number "38a" to decimal
base_11_number = "38a"
decimal_value = 0
base_11_digits = {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'a': 10}

# Calculate the decimal value
for i, digit in enumerate(reversed(base_11_number)):
    decimal_value += base_11_digits[digit] * (11 ** i)

# Step 2: Convert the decimal number to base-8
base_8_number = ""
while decimal_value > 0:
    remainder = decimal_value % 8
    base_8_number = str(remainder) + base_8_number
    decimal_value //= 8

print(base_8_number)