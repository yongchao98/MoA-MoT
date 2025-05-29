# Step 1: Convert base-11 number "20a" to decimal
base_11_number = "20a"
base_11 = 11

# Convert base-11 to decimal
decimal_value = 0
for index, digit in enumerate(reversed(base_11_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_11 ** index)

# Step 2: Convert decimal to base-14
base_14 = 14
base_14_number = []

while decimal_value > 0:
    remainder = decimal_value % base_14
    if remainder >= 10:
        base_14_number.append(chr(remainder - 10 + ord('a')))
    else:
        base_14_number.append(str(remainder))
    decimal_value //= base_14

# The base-14 number is the reverse of the remainders collected
base_14_number.reverse()
base_14_result = ''.join(base_14_number)

print(base_14_result)