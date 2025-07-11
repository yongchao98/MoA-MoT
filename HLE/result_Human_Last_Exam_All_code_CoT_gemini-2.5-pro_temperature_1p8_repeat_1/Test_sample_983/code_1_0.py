def solve_sequence():
    """
    Calculates the next number in a sequence based on a derived rule.
    The rule is: next_number = current_number + sum_of_digits + product_of_digits.
    """
    current_number = 2352

    # Convert the number to a string to iterate through its digits
    s_num = str(current_number)
    
    # Calculate the sum of the digits
    sum_of_digits = 0
    for digit in s_num:
        sum_of_digits += int(digit)
        
    # Calculate the product of the digits
    product_of_digits = 1
    for digit in s_num:
        product_of_digits *= int(digit)
        
    # Calculate the number to add
    added_value = sum_of_digits + product_of_digits
    
    # Calculate the next number in the sequence
    next_number = current_number + added_value
    
    # Print the explanation and the final equation
    print(f"The last number in the sequence is {current_number}.")
    print("The rule to find the next number is to add the sum of the digits and the product of the digits to the current number.")
    print(f"Sum of digits of {current_number} is: " + " + ".join(s_num) + f" = {sum_of_digits}")
    print(f"Product of digits of {current_number} is: " + " * ".join(s_num) + f" = {product_of_digits}")
    print(f"The value to add is {sum_of_digits} + {product_of_digits} = {added_value}.")
    print(f"So, the next number is {current_number} + {added_value} = {next_number}.")
    print(f"\nFinal Equation: {current_number} + ({' + '.join(s_num)}) + ({' * '.join(s_num)}) = {next_number}")

solve_sequence()
<<<2424>>>