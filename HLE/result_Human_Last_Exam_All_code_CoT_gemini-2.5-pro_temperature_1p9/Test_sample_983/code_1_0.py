def get_next_number_in_sequence(n):
    """
    Calculates the next number in the sequence based on the rule:
    next = current + sum_of_digits(current) + product_of_digits(current).
    """
    s = str(n)
    digits = [int(digit) for digit in s]
    
    # Calculate sum of digits
    sum_of_digits = sum(digits)
    
    # Calculate product of digits
    product_of_digits = 1
    for digit in digits:
        product_of_digits *= digit
        
    # Calculate the next number
    next_num = n + sum_of_digits + product_of_digits
    
    # Print the equation
    digit_sum_str = " + ".join(str(d) for d in digits)
    digit_prod_str = " * ".join(str(d) for d in digits)
    print(f"The last number is {n}")
    print(f"Sum of its digits is: {digit_sum_str} = {sum_of_digits}")
    print(f"Product of its digits is: {digit_prod_str} = {product_of_digits}")
    print(f"The next number is: {n} + {sum_of_digits} + {product_of_digits} = {next_num}")
    
    return next_num

# The last number in the provided sequence
last_number = 2352

# Get the next number and print the explanation
get_next_number_in_sequence(last_number)
