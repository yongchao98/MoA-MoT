# Step 1: Convert base-8 number 1024 to base-10
base_8_number = "1024"
base_10_number = int(base_8_number, 8)

# Step 2: Convert base-10 number to base-6
def convert_to_base_6(n):
    if n == 0:
        return "0"
    digits = []
    while n:
        digits.append(int(n % 6))
        n //= 6
    return ''.join(str(x) for x in digits[::-1])

base_6_number = convert_to_base_6(base_10_number)
print(base_6_number)