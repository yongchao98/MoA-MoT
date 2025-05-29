# Step 1: Convert base-6 to base-10
base_6_number = '4125'
base_10_number = sum(int(digit) * (6 ** idx) for idx, digit in enumerate(reversed(base_6_number)))

# Step 2: Convert base-10 to base-5
base_5_number = ''
while base_10_number > 0:
    base_5_number = str(base_10_number % 5) + base_5_number
    base_10_number //= 5

print(base_5_number)