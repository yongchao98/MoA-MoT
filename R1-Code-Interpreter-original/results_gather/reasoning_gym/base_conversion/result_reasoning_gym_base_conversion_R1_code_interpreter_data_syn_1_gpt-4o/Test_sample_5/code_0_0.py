def base11_to_decimal(base11_num):
    # Convert base-11 number to decimal
    decimal_num = 0
    base11_str = str(base11_num)
    length = len(base11_str)
    for i, digit in enumerate(base11_str):
        decimal_num += int(digit) * (11 ** (length - i - 1))
    return decimal_num

def decimal_to_base3(decimal_num):
    # Convert decimal number to base-3
    if decimal_num == 0:
        return "0"
    base3_num = ""
    while decimal_num > 0:
        remainder = decimal_num % 3
        base3_num = str(remainder) + base3_num
        decimal_num = decimal_num // 3
    return base3_num

# Convert base-11 number 585 to base-3
base11_num = 585
decimal_num = base11_to_decimal(base11_num)
base3_num = decimal_to_base3(decimal_num)
print(base3_num)