# Step 1: Convert base-14 to decimal
base_14_number = "a8"
base_14 = 14

# Convert each digit to decimal
decimal_value = 0
for i, digit in enumerate(reversed(base_14_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_14 ** i)

# Step 2: Convert decimal to base-9
base_9 = 9
base_9_number = ""

while decimal_value > 0:
    remainder = decimal_value % base_9
    base_9_number = str(remainder) + base_9_number
    decimal_value //= base_9

print(base_9_number)