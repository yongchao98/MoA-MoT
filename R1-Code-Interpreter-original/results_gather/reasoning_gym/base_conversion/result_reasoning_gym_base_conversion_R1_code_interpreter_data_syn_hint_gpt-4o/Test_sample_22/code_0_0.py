# Step 1: Convert base-4 number to base-10
base_4_number = '32213'
base_10_number = sum(int(digit) * (4 ** idx) for idx, digit in enumerate(reversed(base_4_number)))

# Step 2: Convert base-10 number to base-9
base_9_number = ''
while base_10_number > 0:
    remainder = base_10_number % 9
    base_9_number = str(remainder) + base_9_number
    base_10_number //= 9

print(base_9_number)