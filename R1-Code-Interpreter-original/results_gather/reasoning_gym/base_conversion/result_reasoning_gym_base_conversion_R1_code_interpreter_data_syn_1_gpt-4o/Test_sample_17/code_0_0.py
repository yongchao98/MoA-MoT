# Step 1: Convert base-7 number 125 to base-10
base_7_number = "125"
base_10_number = 0
base_7_length = len(base_7_number)

for i in range(base_7_length):
    digit = int(base_7_number[i])
    power = base_7_length - i - 1
    base_10_number += digit * (7 ** power)

# Step 2: Convert base-10 number to base-8
base_8_number = ""
while base_10_number > 0:
    remainder = base_10_number % 8
    base_8_number = str(remainder) + base_8_number
    base_10_number //= 8

print(base_8_number)