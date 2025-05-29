# Step 1: Convert base-3 number to base-10
base_3_number = "1020200"
base_10_number = 0

# Convert base-3 to base-10
for index, digit in enumerate(reversed(base_3_number)):
    base_10_number += int(digit) * (3 ** index)

# Step 2: Convert base-10 number to base-8
base_8_number = ""
while base_10_number > 0:
    remainder = base_10_number % 8
    base_8_number = str(remainder) + base_8_number
    base_10_number //= 8

# Print the result
print(base_8_number)