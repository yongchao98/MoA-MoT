# Step 1: Convert base-8 number to base-10
base_8_number = '307'
base_10_number = int(base_8_number, 8)

# Step 2: Convert base-10 number to base-14
def base_10_to_base_14(n):
    if n == 0:
        return '0'
    digits = []
    while n:
        remainder = n % 14
        if remainder >= 10:
            # Convert remainder to corresponding letter
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        n //= 14
    return ''.join(reversed(digits))

base_14_number = base_10_to_base_14(base_10_number)
print(base_14_number)