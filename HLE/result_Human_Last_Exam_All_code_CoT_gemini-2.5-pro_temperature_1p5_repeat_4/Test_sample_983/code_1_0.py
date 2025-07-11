def solve():
    """
    Calculates the next number in a sequence based on a specific rule.
    The rule is: a_{n+1} = a_n + sum_of_digits(a_n) + product_of_digits(a_n).
    """
    current_number = 2352
    
    # Convert the number to a string to iterate through its digits.
    s_number = str(current_number)
    
    # Calculate the sum of digits.
    sum_of_digits = 0
    for digit in s_number:
        sum_of_digits += int(digit)
        
    # Calculate the product of digits.
    product_of_digits = 1
    for digit in s_number:
        product_of_digits *= int(digit)
        
    # Calculate the next number.
    next_number = current_number + sum_of_digits + product_of_digits
    
    # Print the equation and the final answer.
    print(f"The last number is {current_number}.")
    print(f"The sum of its digits is: " + " + ".join(s_number) + f" = {sum_of_digits}")
    print(f"The product of its digits is: " + " * ".join(s_number) + f" = {product_of_digits}")
    print(f"The next number is calculated as: {current_number} + {sum_of_digits} + {product_of_digits} = {next_number}")

solve()
<<<2424>>>