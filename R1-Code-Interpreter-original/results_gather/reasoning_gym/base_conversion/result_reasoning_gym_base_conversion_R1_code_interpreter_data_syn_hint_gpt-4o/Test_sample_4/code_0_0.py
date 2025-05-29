# Step 1: Convert binary to decimal
binary_number = '101010000'
decimal_number = int(binary_number, 2)

# Step 2: Convert decimal to base-14
def decimal_to_base_n(decimal, base):
    if decimal == 0:
        return '0'
    digits = []
    while decimal:
        remainder = decimal % base
        if remainder >= 10:
            # Convert remainder to corresponding letter for bases > 10
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        decimal //= base
    return ''.join(reversed(digits))

base_14_number = decimal_to_base_n(decimal_number, 14)
print(base_14_number)