# Step 1: Convert base-9 number 1033 to base-10
base_9_number = "1033"
base_10_number = 0
base_9_length = len(base_9_number)

for i in range(base_9_length):
    digit = int(base_9_number[i])
    power = base_9_length - i - 1
    base_10_number += digit * (9 ** power)

# Step 2: Convert base-10 number to base-13
base_13_number = ""
while base_10_number > 0:
    remainder = base_10_number % 13
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_13_number = chr(remainder - 10 + ord('a')) + base_13_number
    else:
        base_13_number = str(remainder) + base_13_number
    base_10_number //= 13

print(base_13_number)