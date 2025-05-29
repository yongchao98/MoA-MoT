# Step 1: Convert base-9 number 555 to base-10
base_9_number = "555"
base_10_number = sum(int(digit) * (9 ** idx) for idx, digit in enumerate(reversed(base_9_number)))

# Step 2: Convert base-10 number to base-6
base_6_number = ""
while base_10_number > 0:
    remainder = base_10_number % 6
    base_6_number = str(remainder) + base_6_number
    base_10_number //= 6

print(base_6_number)