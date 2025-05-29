# Step 1: Convert hexadecimal '2d' to decimal
hex_number = '2d'
decimal_number = int(hex_number, 16)

# Step 2: Convert decimal to base-11
def decimal_to_base_n(decimal, base):
    if decimal == 0:
        return '0'
    digits = []
    while decimal:
        remainder = decimal % base
        if remainder >= 10:
            # Convert remainder to corresponding letter for base > 10
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        decimal //= base
    return ''.join(reversed(digits))

base_11_number = decimal_to_base_n(decimal_number, 11)
print(base_11_number)