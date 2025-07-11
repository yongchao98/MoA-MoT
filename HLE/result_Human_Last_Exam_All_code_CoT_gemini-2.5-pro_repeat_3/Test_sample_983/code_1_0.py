def solve_sequence():
    """
    Calculates the next number in the sequence by adding the sum of the digits
    of the last number to itself.
    """
    last_number = 2352
    
    # Calculate the sum of the digits of the last number
    sum_of_digits = sum(int(digit) for digit in str(last_number))
    
    # Calculate the next number
    next_number = last_number + sum_of_digits
    
    # Print the equation as requested
    print(f"{last_number} + {sum_of_digits} = {next_number}")

solve_sequence()