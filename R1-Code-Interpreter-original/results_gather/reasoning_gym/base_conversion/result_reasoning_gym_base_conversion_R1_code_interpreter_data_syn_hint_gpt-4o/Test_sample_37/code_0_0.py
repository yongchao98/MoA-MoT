# Step 1: Convert binary to decimal
binary_number = '1010111100'
decimal_number = int(binary_number, 2)

# Step 2: Convert decimal to base 3
def decimal_to_base3(n):
    if n == 0:
        return '0'
    base3 = ''
    while n > 0:
        base3 = str(n % 3) + base3
        n = n // 3
    return base3

# Convert the decimal number to base 3
base3_number = decimal_to_base3(decimal_number)

print(base3_number)