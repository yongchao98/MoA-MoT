# Step 1: Convert base-4 to base-10
base_4_number = '20233'
base_10_number = sum(int(digit) * (4 ** idx) for idx, digit in enumerate(reversed(base_4_number)))

# Step 2: Convert base-10 to base-3
base_3_number = ''
while base_10_number > 0:
    base_3_number = str(base_10_number % 3) + base_3_number
    base_10_number //= 3

print(base_3_number)