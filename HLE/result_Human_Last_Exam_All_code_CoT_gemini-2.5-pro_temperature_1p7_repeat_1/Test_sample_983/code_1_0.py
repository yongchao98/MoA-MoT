def solve_sequence():
    """
    This function calculates the next number in the sequence based on the rule:
    next_number = current_number + sum_of_digits(current_number).
    """
    
    # The last number in the given sequence
    last_number = 2352
    
    # Calculate the sum of the digits of the last number
    current_number_str = str(last_number)
    sum_of_digits = 0
    digit_sum_str_parts = []
    for digit in current_number_str:
        sum_of_digits += int(digit)
        digit_sum_str_parts.append(digit)

    digit_sum_str = " + ".join(digit_sum_str_parts)

    # Calculate the next number in the sequence
    next_number = last_number + sum_of_digits
    
    # Print the equation explaining the result
    print(f"The last number is {last_number}.")
    print(f"The sum of its digits is {digit_sum_str} = {sum_of_digits}.")
    print(f"The rule is to add the sum of the digits to the number itself.")
    print(f"So, the next number is {last_number} + {sum_of_digits} = {next_number}.")

solve_sequence()