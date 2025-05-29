# Step 1: Convert base-6 number 200 to base-10
base_6_number = "200"
base_10_number = int(base_6_number, 6)

# Step 2: Convert base-10 number to base-8
base_8_number = ""
while base_10_number > 0:
    remainder = base_10_number % 8
    base_8_number = str(remainder) + base_8_number
    base_10_number //= 8

print(base_8_number)