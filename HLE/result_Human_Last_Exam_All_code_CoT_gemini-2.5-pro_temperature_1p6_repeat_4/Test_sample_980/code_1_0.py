def find_next_number():
    """
    Calculates the next number based on the sum of the squares of the digits
    of the last number in the provided sequence.
    """
    last_number = 124
    
    # Convert the number to a string to iterate through its digits.
    s_number = str(last_number)
    
    # Create lists to hold the parts of the equation.
    equation_parts = []
    digit_values = []
    
    for digit_char in s_number:
        digit = int(digit_char)
        equation_parts.append(f"{digit}^2")
        digit_values.append(digit)
        
    # Calculate the final result.
    result = sum(d**2 for d in digit_values)
    
    # Format the equation string.
    equation_string = " + ".join(equation_parts)
    
    # Print the full equation as requested.
    print(f"The final equation is: {equation_string} = {result}")

find_next_number()