def get_next_number(n):
    """
    Calculates the next number in the sequence based on the rule:
    n + (product of non-zero digits of n) / (first digit of n)
    """
    s = str(n)
    
    # Find the first digit
    first_digit = int(s[0])
    
    # Calculate the product of non-zero digits
    product_nz = 1
    for digit in s:
        if digit != '0':
            product_nz *= int(digit)
            
    # Calculate the value to add
    value_to_add = product_nz // first_digit
    
    # Calculate the next number
    next_n = n + value_to_add
    
    # Print the equation with all the numbers
    print(f"The last number is {n}.")
    print(f"The digits of {n} are {', '.join(list(s))}.")
    print(f"The first digit is {first_digit}.")
    print(f"The product of its non-zero digits is {' * '.join([d for d in s if d != '0'])} = {product_nz}.")
    print(f"The value to add is the product divided by the first digit: {product_nz} / {first_digit} = {value_to_add}.")
    print(f"The final equation is: {n} + {value_to_add} = {next_n}")
    print(f"The next number in the sequence is {next_n}.")

# The last number in the sequence provided by the user
last_number = 2352
get_next_number(last_number)
