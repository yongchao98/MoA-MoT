def solve():
    """
    This function calculates the next number in the sequence based on the rule:
    next_number = current_number + 6 * (sum of the two middle digits).
    """
    last_number = 2352
    
    # Extract the two middle digits (hundreds and tens place)
    # For a number n=abcd, hundreds_digit is b, tens_digit is c
    hundreds_digit = (last_number // 100) % 10
    tens_digit = (last_number // 10) % 10
    
    # Calculate the increment based on the rule
    increment = 6 * (hundreds_digit + tens_digit)
    
    # Calculate the next number in the sequence
    next_number = last_number + increment
    
    # Print the equation
    print(f"The last number is: {last_number}")
    print(f"The rule is: current_number + 6 * (hundreds_digit + tens_digit)")
    print(f"Calculation: {last_number} + 6 * ({hundreds_digit} + {tens_digit}) = {last_number} + 6 * {hundreds_digit + tens_digit} = {last_number} + {increment} = {next_number}")
    print(f"The next number in the sequence is: {next_number}")

solve()
