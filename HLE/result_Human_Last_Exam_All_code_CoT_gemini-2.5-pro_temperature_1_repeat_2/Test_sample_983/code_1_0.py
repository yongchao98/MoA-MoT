def solve_sequence():
    """
    Calculates the next number in a sequence based on a derived rule.
    """
    num = 2352
    
    # Convert the number to a list of its digits
    digits = [int(d) for d in str(num)]
    
    # Calculate the sum of the digits
    sum_of_digits = sum(digits)
    
    # Calculate the product of the digits
    product_of_digits = 1
    for digit in digits:
        product_of_digits *= digit
        
    # Calculate the next number in the sequence
    next_num = num + sum_of_digits + product_of_digits
    
    # Create string representations for printing the equation
    s_o_d_str = " + ".join(map(str, digits))
    p_o_d_str = " * ".join(map(str, digits))
    
    print(f"To find the next number, we apply the rule: number + sum_of_digits + product_of_digits")
    print(f"Last Number: {num}")
    print(f"Sum of Digits: {s_o_d_str} = {sum_of_digits}")
    print(f"Product of Digits: {p_o_d_str} = {product_of_digits}")
    print(f"The final equation is: {num} + {sum_of_digits} + {product_of_digits} = {next_num}")
    print(f"The next number in the sequence is: {next_num}")

solve_sequence()