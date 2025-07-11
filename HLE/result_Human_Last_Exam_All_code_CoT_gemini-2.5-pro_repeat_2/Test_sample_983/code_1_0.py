def sum_digits(n):
    """Calculates the sum of the digits of a number."""
    s = 0
    for digit in str(n):
        s += int(digit)
    return s

def product_digits(n):
    """Calculates the product of the digits of a number."""
    if '0' in str(n):
        return 0
    p = 1
    for digit in str(n):
        p *= int(digit)
    return p

# The given sequence
sequence = [1014, 1020, 1032, 1050, 1104, 1224, 1230, 1290, 1410, 1440, 1452, 1482, 1500, 1512, 1560, 1584, 1602, 1734, 1812, 1890, 1902, 1932, 1974, 2004, 2040, 2064, 2070, 2130, 2142, 2274, 2340, 2352]

# The last number in the sequence
last_number = sequence[-1]

# Calculate the sum and product of its digits
s_digits = sum_digits(last_number)
p_digits = product_digits(last_number)

# Apply the rule: next_number = current_number + sum_of_digits + product_of_digits
next_number = last_number + s_digits + p_digits

# Output the equation and the result
print(f"The last number is {last_number}.")
print(f"The sum of its digits is {' + '.join(str(last_number))} = {s_digits}.")
print(f"The product of its digits is {' * '.join(str(last_number))} = {p_digits}.")
print(f"The rule is: next_number = current_number + sum_of_digits + product_of_digits.")
print(f"So, the next number is {last_number} + {s_digits} + {p_digits} = {next_number}")
