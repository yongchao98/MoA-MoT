def sum_digits(n):
    """Calculates the sum of the digits of a number."""
    s = 0
    for digit in str(n):
        s += int(digit)
    return s

def product_nonzero_digits(n):
    """Calculates the product of the non-zero digits of a number."""
    p = 1
    for digit in str(n):
        if int(digit) != 0:
            p *= int(digit)
    return p

def find_next_number(n):
    """
    Calculates the next number in the sequence based on the rule:
    n_next = n + (sum_digits(n) * product_nonzero_digits(n))
    """
    s = sum_digits(n)
    p = product_nonzero_digits(n)
    next_n = n + s * p
    print(f"The last number is {n}.")
    print(f"The sum of its digits is: {' + '.join(str(n))} = {s}")
    print(f"The product of its non-zero digits is: {' * '.join(c for c in str(n) if c != '0')} = {p}")
    print(f"The number to add is the product of these two: {s} * {p} = {s*p}")
    print(f"The next number in the sequence is: {n} + {s*p} = {next_n}")
    return next_n

# The last number in the provided sequence
last_number = 2352
final_answer = find_next_number(last_number)
