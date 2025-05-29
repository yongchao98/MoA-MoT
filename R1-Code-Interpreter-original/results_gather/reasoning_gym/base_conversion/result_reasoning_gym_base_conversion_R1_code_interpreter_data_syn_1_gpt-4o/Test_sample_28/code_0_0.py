# Convert base-6 number 1243 to base-10
base_6_number = '1243'
base_10_number = sum(int(digit) * (6 ** idx) for idx, digit in enumerate(reversed(base_6_number)))

# Convert base-10 number to base-11
base_11_number = ''
while base_10_number > 0:
    remainder = base_10_number % 11
    if remainder >= 10:
        # Convert remainder to corresponding letter for bases > 10
        base_11_number = chr(remainder - 10 + ord('a')) + base_11_number
    else:
        base_11_number = str(remainder) + base_11_number
    base_10_number //= 11

print(base_11_number)