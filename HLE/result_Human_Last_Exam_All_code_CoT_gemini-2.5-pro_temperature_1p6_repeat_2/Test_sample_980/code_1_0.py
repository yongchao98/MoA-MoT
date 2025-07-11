def solve():
    """
    Calculates the next number in the sequence by applying a transformation
    to the last number. The transformation is the sum of the squares of the digits.
    """
    last_number_str = "124"
    
    # Get the digits of the last number
    digits = [int(d) for d in last_number_str]
    
    # Calculate the squares of the digits
    squares = [d * d for d in digits]
    
    # Sum the squares to get the result
    result = sum(squares)
    
    # Build the equation string part by part to meet the output requirement
    equation_parts = []
    for d in digits:
        equation_parts.append(f"{d}*{d}")
    
    equation_str = " + ".join(equation_parts)
    
    # Print the final equation with the result
    print(f"{equation_str} = {result}")

solve()