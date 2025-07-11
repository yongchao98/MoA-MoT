def solve():
    """
    This function calculates the next number in the sequence based on a discovered pattern.
    """
    # The last number in the provided sequence.
    current_number = 2352

    # Calculate the sum of the digits of the current number.
    s = str(current_number)
    sum_of_digits = sum(int(digit) for digit in s)

    # The pattern suggests the next number is the current number plus the sum of its digits.
    next_number = current_number + sum_of_digits

    # Print the equation showing how the next number is calculated.
    # The list comprehension creates a list of strings from the digits, then join adds '+' between them.
    digit_sum_str = ' + '.join(list(s))
    print(f"The next number is found by adding the sum of the digits of the last number to itself.")
    print(f"The last number is {current_number}.")
    print(f"The sum of its digits is {digit_sum_str} = {sum_of_digits}.")
    print(f"So, the final equation is: {current_number} + {sum_of_digits} = {next_number}")

solve()
<<<2364>>>