# Step 1: Convert base-12 number "27" to decimal
base_12_number = "27"
decimal_value = int(base_12_number, 12)

# Step 2: Convert the decimal value to base-15
def decimal_to_base_n(decimal_number, base):
    if decimal_number == 0:
        return "0"
    digits = []
    while decimal_number:
        remainder = decimal_number % base
        if remainder >= 10:
            # Convert remainder to corresponding letter for bases > 10
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        decimal_number //= base
    return ''.join(reversed(digits))

base_15_number = decimal_to_base_n(decimal_value, 15)
print(base_15_number)