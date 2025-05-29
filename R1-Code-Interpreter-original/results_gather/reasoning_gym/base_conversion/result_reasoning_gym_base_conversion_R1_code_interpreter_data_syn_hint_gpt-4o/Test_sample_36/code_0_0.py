# Step 1: Convert base-9 number 1033 to base-10
base_9_number = "1033"
base_10_number = sum(int(digit) * (9 ** idx) for idx, digit in enumerate(reversed(base_9_number)))

# Step 2: Convert base-10 number to base-13
base_13_number = ""
base = 13
while base_10_number > 0:
    remainder = base_10_number % base
    if remainder >= 10:
        # Convert remainder to corresponding letter for bases > 10
        base_13_number = chr(remainder - 10 + ord('a')) + base_13_number
    else:
        base_13_number = str(remainder) + base_13_number
    base_10_number //= base

print(base_13_number)