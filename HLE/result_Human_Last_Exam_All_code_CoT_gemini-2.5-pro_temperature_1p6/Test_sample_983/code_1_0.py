def solve():
    """
    This function calculates the next number in a sequence where the increment
    is derived from the digits of the current number.
    """
    # The last number in the given sequence
    last_number = 2352

    # Step 1: Calculate the sum of the digits of the last number.
    # The digits are 2, 3, 5, 2.
    sum_of_digits = 2 + 3 + 5 + 2

    # Step 2: For this specific step in the sequence, we determine the
    # increment is half of the sum of the digits.
    increment = sum_of_digits // 2

    # Step 3: Calculate the next number by adding the increment to the last number.
    next_number = last_number + increment

    # Output the final equation showing each number.
    print(f"The final calculation is: {last_number} + {increment} = {next_number}")


solve()