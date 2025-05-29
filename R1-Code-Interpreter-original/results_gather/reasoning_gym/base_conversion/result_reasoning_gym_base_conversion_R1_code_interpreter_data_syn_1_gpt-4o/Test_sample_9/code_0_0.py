# Step 1: Convert base-7 to base-10
base_7_number = "1065"
decimal_number = 0
base_7_length = len(base_7_number)

for i in range(base_7_length):
    digit = int(base_7_number[i])
    power = base_7_length - i - 1
    decimal_number += digit * (7 ** power)

# Step 2: Convert base-10 to base-5
base_5_number = ""
while decimal_number > 0:
    remainder = decimal_number % 5
    base_5_number = str(remainder) + base_5_number
    decimal_number //= 5

print(base_5_number)