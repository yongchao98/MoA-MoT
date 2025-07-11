def solve():
    """
    This function calculates the next number in the sequence.
    The rule is: next_number = current_number + product_of_nonzero_digits(current_number).
    """
    last_num = 2352

    # Find the non-zero digits
    s_num = str(last_num)
    product = 1
    nonzero_digits = []
    for digit in s_num:
        if digit != '0':
            d = int(digit)
            product *= d
            nonzero_digits.append(str(d))

    next_num = last_num + product

    # Format the equation string
    # e.g., "2352 + 2 * 3 * 5 * 2 = 2412"
    equation_str = f"{last_num} + {' * '.join(nonzero_digits)} = {next_num}"
    
    print(equation_str)


solve()