# Step 1: Convert base-5 number 3112 to base-10
base_5_number = '3112'
base_10_number = sum(int(digit) * (5 ** idx) for idx, digit in enumerate(reversed(base_5_number)))

# Step 2: Convert base-10 number to base-7
base_7_number = ''
while base_10_number > 0:
    remainder = base_10_number % 7
    base_7_number = str(remainder) + base_7_number
    base_10_number //= 7

print(base_7_number)