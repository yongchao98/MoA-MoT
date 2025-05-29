# Step 1: Convert base-7 to base-10
base_7_number = "1431"
base_10_number = sum(int(digit) * (7 ** idx) for idx, digit in enumerate(reversed(base_7_number)))

# Step 2: Convert base-10 to base-4
base_4_number = ""
while base_10_number > 0:
    base_4_number = str(base_10_number % 4) + base_4_number
    base_10_number //= 4

print(base_4_number)