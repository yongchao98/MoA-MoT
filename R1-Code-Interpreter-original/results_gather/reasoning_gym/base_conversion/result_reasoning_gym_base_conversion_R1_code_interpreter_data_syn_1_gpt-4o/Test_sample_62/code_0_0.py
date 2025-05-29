# Step 1: Convert base-4 number 20233 to base-10
base_4_number = "20233"
base_10_number = 0
base_4_length = len(base_4_number)

for i in range(base_4_length):
    digit = int(base_4_number[i])
    power = base_4_length - i - 1
    base_10_number += digit * (4 ** power)

# Step 2: Convert base-10 number to base-3
base_3_number = ""
while base_10_number > 0:
    remainder = base_10_number % 3
    base_3_number = str(remainder) + base_3_number
    base_10_number //= 3

# If the number is 0, we should return "0"
if base_3_number == "":
    base_3_number = "0"

print(base_3_number)