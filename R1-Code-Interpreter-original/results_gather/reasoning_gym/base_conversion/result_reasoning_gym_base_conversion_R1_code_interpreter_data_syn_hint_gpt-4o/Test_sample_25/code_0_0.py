# Step 1: Convert base-7 number 61 to decimal
base_7_number = "61"
decimal_number = int(base_7_number, 7)

# Step 2: Convert decimal number to base-11
def decimal_to_base_n(decimal, base):
    if decimal == 0:
        return "0"
    digits = []
    while decimal:
        digits.append(int(decimal % base))
        decimal //= base
    # Convert digits to appropriate base-11 representation
    base_n_number = ''.join(str(d) if d < 10 else chr(d - 10 + ord('a')) for d in reversed(digits))
    return base_n_number

base_11_number = decimal_to_base_n(decimal_number, 11)
print(base_11_number)