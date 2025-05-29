# Convert base-5 number 4334 to base-10
base_5_number = '4334'
base_5 = 5

# Convert to base-10
base_10_number = sum(int(digit) * (base_5 ** idx) for idx, digit in enumerate(reversed(base_5_number)))

# Convert base-10 number to base-13
base_13 = 13
base_13_number = ''
while base_10_number > 0:
    remainder = base_10_number % base_13
    if remainder >= 10:
        # Convert remainder to corresponding letter
        base_13_number = chr(remainder - 10 + ord('a')) + base_13_number
    else:
        base_13_number = str(remainder) + base_13_number
    base_10_number //= base_13

print(base_13_number)