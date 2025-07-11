def find_next_number():
    """
    Calculates the next number in a sequence based on the rule:
    next_number = current_number + sum_of_its_digits.
    """
    # The last number given in the sequence.
    current_number = 124

    # Convert the number to a string to iterate through its digits.
    s_number = str(current_number)

    # Calculate the sum of the digits.
    sum_of_digits = 0
    for digit_char in s_number:
        sum_of_digits += int(digit_char)

    # Calculate the next number.
    next_number = current_number + sum_of_digits

    # Create the equation string, showing each number involved.
    # The digits are joined by ' + ' for a clear representation of the sum.
    digit_sum_str = " + ".join(s_number)
    
    # Print the final equation.
    print(f"{current_number} + {digit_sum_str} = {next_number}")

find_next_number()