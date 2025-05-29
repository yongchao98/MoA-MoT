# Step 1: Convert base-15 to decimal
base_15_number = "22c"
base_15 = 15

# Convert each digit to its decimal equivalent
decimal_value = 0
for i, digit in enumerate(reversed(base_15_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_15 ** i)

# Step 2: Convert decimal to base-12
base_12 = 12
base_12_number = []

while decimal_value > 0:
    remainder = decimal_value % base_12
    if remainder >= 10:
        base_12_number.append(chr(remainder - 10 + ord('a')))
    else:
        base_12_number.append(str(remainder))
    decimal_value //= base_12

# The base-12 number is the reverse of the remainders collected
base_12_number.reverse()
base_12_result = ''.join(base_12_number)

print(base_12_result)