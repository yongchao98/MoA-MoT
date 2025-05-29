# Step 1: Convert base-14 number '41c' to decimal
base_14_number = '41c'
base_14 = 14

# Convert each digit to its decimal equivalent
decimal_value = 0
for index, digit in enumerate(reversed(base_14_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_14 ** index)

# Step 2: Convert the decimal number to binary
binary_value = ''
while decimal_value > 0:
    binary_value = str(decimal_value % 2) + binary_value
    decimal_value //= 2

print(binary_value)