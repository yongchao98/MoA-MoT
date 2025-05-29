# Step 1: Convert base-13 number 'b6' to decimal
base_13_number = 'b6'
# Calculate the decimal value manually
decimal_value = 11 * 13 + 6

# Step 2: Convert the decimal value to base-5
def decimal_to_base_n(decimal_number, base):
    if decimal_number == 0:
        return '0'
    digits = []
    while decimal_number:
        digits.append(int(decimal_number % base))
        decimal_number //= base
    return ''.join(str(x) for x in digits[::-1])

# Convert the decimal value to base-5
base_5_value = decimal_to_base_n(decimal_value, 5)

print(base_5_value)