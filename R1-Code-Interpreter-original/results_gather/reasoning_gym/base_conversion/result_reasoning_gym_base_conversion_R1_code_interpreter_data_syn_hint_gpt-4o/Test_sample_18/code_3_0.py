# Step 1: Convert base-15 to decimal
base_15_number = "23e"
base_15 = 15

# Convert each digit to its decimal equivalent
decimal_value = 0
for i, digit in enumerate(reversed(base_15_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_15 ** i)

# Step 2: Convert decimal to base-6
base_6 = 6
base_6_number = []

while decimal_value > 0:
    remainder = decimal_value % base_6
    base_6_number.append(str(remainder))
    decimal_value //= base_6

# The base-6 number is the remainders read in reverse order
base_6_number.reverse()
base_6_result = ''.join(base_6_number)

print(base_6_result)