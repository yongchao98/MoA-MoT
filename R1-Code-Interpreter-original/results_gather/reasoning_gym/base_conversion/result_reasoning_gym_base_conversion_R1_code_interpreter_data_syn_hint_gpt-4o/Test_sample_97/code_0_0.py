# Step 1: Convert base-15 number "2a0" to decimal
base_15_number = "2a0"
decimal_value = 0
base_15_digits = {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'a': 10, 'b': 11, 'c': 12, 'd': 13, 'e': 14}

for i, digit in enumerate(reversed(base_15_number)):
    decimal_value += base_15_digits[digit] * (15 ** i)

# Step 2: Convert decimal number to base-14
base_14_digits = "0123456789abcd"
base_14_number = ""

while decimal_value > 0:
    remainder = decimal_value % 14
    base_14_number = base_14_digits[remainder] + base_14_number
    decimal_value //= 14

print(base_14_number)