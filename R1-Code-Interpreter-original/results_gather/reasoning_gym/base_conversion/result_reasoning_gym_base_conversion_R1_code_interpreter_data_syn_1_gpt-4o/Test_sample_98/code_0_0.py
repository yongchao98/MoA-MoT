# Step 1: Convert base-8 number 1024 to base-10
base_8_number = "1024"
base_10_number = 0
base_8_length = len(base_8_number)

for i in range(base_8_length):
    digit = int(base_8_number[i])
    power = base_8_length - i - 1
    base_10_number += digit * (8 ** power)

# Step 2: Convert base-10 number to base-6
base_6_number = ""
current_number = base_10_number

while current_number > 0:
    remainder = current_number % 6
    base_6_number = str(remainder) + base_6_number
    current_number //= 6

print(base_6_number)