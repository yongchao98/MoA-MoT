def solve_sequence():
    """
    Calculates the next number in the sequence based on the observed pattern at the end.
    """
    last_number = 2352

    # Convert the number to a string to iterate through its digits
    s = str(last_number)
    
    # Calculate the sum of the digits
    sum_of_digits = sum(int(digit) for digit in s)
    
    # Calculate the next number
    next_number = last_number + sum_of_digits
    
    # Print the equation step by step
    digit_sum_str = " + ".join(s)
    print(f"The last number is {last_number}.")
    print(f"The sum of its digits is: {digit_sum_str} = {sum_of_digits}.")
    print(f"The equation for the next number is: {last_number} + {sum_of_digits} = {next_number}")
    print(f"The next number in the sequence is {next_number}.")

solve_sequence()
<<<2364>>>