def solve():
    """
    Calculates the next number in the sequence based on the rule:
    next_num = current_num + (sum of digits)^2 - (sum of squares of digits).
    """
    current_num = 2352
    
    # Extract digits from the number
    s_num = str(current_num)
    digits = [int(d) for d in s_num]
    
    # Calculate sum of digits (S)
    sum_of_digits = sum(digits)
    
    # Calculate sum of squares of digits
    sum_of_squares = sum(d**2 for d in digits)
    
    # Calculate the value to be added
    added_value = sum_of_digits**2 - sum_of_squares
    
    # Calculate the next number in the sequence
    next_num = current_num + added_value
    
    print(f"The last number in the sequence is {current_num}.")
    print(f"The digits are: {', '.join(map(str, digits))}")
    print(f"The sum of the digits is: {' + '.join(map(str, digits))} = {sum_of_digits}")
    print(f"The sum of the squares of the digits is: {' + '.join([f'{d}^2' for d in digits])} = {sum_of_squares}")
    print(f"The rule is: next = current + (sum_of_digits^2 - sum_of_squares_of_digits)")
    print(f"The value to add is: {sum_of_digits}^2 - {sum_of_squares} = {added_value}")
    print(f"The final equation is: {current_num} + {added_value} = {next_num}")
    print(f"The next number in the sequence is {next_num}.")

solve()