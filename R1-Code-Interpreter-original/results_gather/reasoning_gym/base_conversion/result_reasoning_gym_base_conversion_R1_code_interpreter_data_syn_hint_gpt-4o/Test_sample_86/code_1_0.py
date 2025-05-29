# Step 1: Convert base-12 to base-10
base_12_number = "508"
base_10_number = 5 * (12**2) + 0 * (12**1) + 8 * (12**0)

# Step 2: Convert base-10 to base-4
def convert_to_base_4(n):
    if n == 0:
        return "0"
    digits = []
    while n:
        digits.append(int(n % 4))
        n //= 4
    return ''.join(str(x) for x in digits[::-1])

base_4_number = convert_to_base_4(base_10_number)
print(base_4_number)