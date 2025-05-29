# Step 1: Convert base-7 number 61 to base-10
base_7_number = "61"
base_10_number = int(base_7_number, 7)

# Step 2: Convert base-10 number to base-11
def convert_to_base_11(n):
    if n == 0:
        return "0"
    digits = []
    while n:
        digits.append(int(n % 11))
        n //= 11
    # Convert digits to appropriate base-11 representation
    base_11_number = ''.join(str(d) if d < 10 else chr(d - 10 + ord('a')) for d in reversed(digits))
    return base_11_number

base_11_number = convert_to_base_11(base_10_number)
print(base_11_number)