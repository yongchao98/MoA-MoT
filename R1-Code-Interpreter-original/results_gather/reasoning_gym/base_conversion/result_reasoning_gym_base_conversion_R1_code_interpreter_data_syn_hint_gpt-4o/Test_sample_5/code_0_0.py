# Step 1: Convert base-11 number 585 to base-10
base_11_number = "585"
base_10_number = int(base_11_number, 11)

# Step 2: Convert base-10 number to base-3
def convert_to_base_3(n):
    if n == 0:
        return "0"
    digits = []
    while n:
        digits.append(int(n % 3))
        n //= 3
    return ''.join(str(x) for x in digits[::-1])

base_3_number = convert_to_base_3(base_10_number)
print(base_3_number)