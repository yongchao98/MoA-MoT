def get_digits(n):
    """Returns a list of digits of a number."""
    return [int(d) for d in str(n)]

def sum_digits(digits):
    """Calculates the sum of a list of digits."""
    return sum(digits)

def product_digits(digits):
    """Calculates the product of a list of digits."""
    p = 1
    for d in digits:
        p *= d
    return p

def calculate_next_number(n):
    """
    Calculates the next number in the sequence based on the rule:
    next_num = current_num + sum_of_digits + product_of_digits.
    """
    digits = get_digits(n)
    sod = sum_digits(digits)
    pod = product_digits(digits)
    increment = sod + pod
    next_num = n + increment
    return next_num, sod, pod, increment

# The last number in the provided sequence
last_number = 2352

# Calculate the next number and the components of the calculation
next_number, sum_val, prod_val, increment_val = calculate_next_number(last_number)

# Print the final equation
print(f"The rule is: next_number = current_number + (sum of digits) + (product of digits)")
print(f"For the number {last_number}:")
print(f"Sum of digits = {' + '.join(str(d) for d in get_digits(last_number))} = {sum_val}")
print(f"Product of digits = {' * '.join(str(d) for d in get_digits(last_number))} = {prod_val}")
print(f"Increment = {sum_val} + {prod_val} = {increment_val}")
print(f"Final Equation: {last_number} + {increment_val} = {next_number}")

print(f"\nThe next number in the sequence is {next_number}")
<<<2424>>>