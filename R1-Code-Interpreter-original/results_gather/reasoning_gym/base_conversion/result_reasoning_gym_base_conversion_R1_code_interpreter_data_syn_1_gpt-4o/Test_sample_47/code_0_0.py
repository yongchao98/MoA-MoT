# Step 1: Convert base-9 number 174 to base-10
base_9_number = "174"
base_10_number = 0
base_9_length = len(base_9_number)

for i in range(base_9_length):
    digit = int(base_9_number[i])
    power = base_9_length - i - 1
    base_10_number += digit * (9 ** power)

# Step 2: Convert base-10 number to base-7
base_7_number = ""
while base_10_number > 0:
    remainder = base_10_number % 7
    base_7_number = str(remainder) + base_7_number
    base_10_number //= 7

print(base_7_number)