def solve():
    """
    This function calculates the next number in a sequence based on a discovered pattern.
    The pattern is that the next number is the current number plus the sum of its digits.
    This pattern is applied to the last number of the given sequence.
    """
    sequence = [
        1014, 1020, 1032, 1050, 1104, 1224, 1230, 1290, 1410, 1440, 1452, 1482, 
        1500, 1512, 1560, 1584, 1602, 1734, 1812, 1890, 1902, 1932, 1974, 2004, 
        2040, 2064, 2070, 2130, 2142, 2274, 2340, 2352
    ]

    # Get the last number in the sequence
    last_number = sequence[-1]

    # Calculate the sum of the digits of the last number
    s = str(last_number)
    sum_of_digits = sum(int(digit) for digit in s)

    # Calculate the next number
    next_number = last_number + sum_of_digits

    # Print the equation
    digit_sum_str = " + ".join(list(s))
    print(f"The last number is {last_number}.")
    print(f"The sum of its digits is: {digit_sum_str} = {sum_of_digits}.")
    print(f"The next number is calculated as: {last_number} + {sum_of_digits} = {next_number}")
    print(f"So, the next number in the sequence is {next_number}.")

solve()
<<<2364>>>