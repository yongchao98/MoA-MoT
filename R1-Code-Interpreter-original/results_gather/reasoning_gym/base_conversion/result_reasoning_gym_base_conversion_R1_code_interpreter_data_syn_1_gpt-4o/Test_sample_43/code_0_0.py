# Step 1: Convert base-11 number "38a" to decimal
base_11_number = "38a"
base_11 = 11

# Convert each character to its decimal equivalent
decimal_value = 0
for i, char in enumerate(reversed(base_11_number)):
    if char.isdigit():
        value = int(char)
    else:
        value = ord(char) - ord('a') + 10  # Convert 'a' to 10, 'b' to 11, etc.
    decimal_value += value * (base_11 ** i)

# Step 2: Convert the decimal number to base-8
base_8 = 8
base_8_number = ""

while decimal_value > 0:
    remainder = decimal_value % base_8
    base_8_number = str(remainder) + base_8_number
    decimal_value //= base_8

print(base_8_number)