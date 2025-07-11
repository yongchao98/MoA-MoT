def solve_sequence():
    """
    Calculates the next number in the sequence by adding the sum of the squares
    of the digits of the last number to itself.
    """
    last_number = 124
    
    # Extract digits
    temp_n = last_number
    digits = []
    while temp_n > 0:
        digits.append(temp_n % 10)
        temp_n //= 10
    digits.reverse()
    
    # Calculate sum of squares of digits
    sum_of_squares = sum(d**2 for d in digits)
    
    # Calculate the next number
    next_number = last_number + sum_of_squares
    
    # Format the equation string as requested
    squared_digits_str = " + ".join([f"{d}**2" for d in digits])
    
    print(f"The last number is {last_number}.")
    print("The rule is: next_number = last_number + (sum of the squares of the digits of last_number).")
    print("Calculation:")
    print(f"{last_number} + ({squared_digits_str}) = {last_number} + {sum_of_squares} = {next_number}")
    print(f"The next number in the sequence is {next_number}.")

solve_sequence()
<<<145>>>