def solve():
    """
    Calculates the next number in the sequence by adding the sum of the digits
    of the last number to itself.
    """
    last_number = 2352
    
    # Calculate the sum of the digits of the last number
    s = str(last_number)
    digits = [int(d) for d in s]
    digit_sum = sum(digits)
    
    # Calculate the next number
    next_number = last_number + digit_sum
    
    # Create the equation string
    # "2352 + 2 + 3 + 5 + 2 = 2364"
    equation_parts = [str(last_number), " + "]
    equation_parts.extend([d + " + " for d in s])
    # remove the last " + "
    equation_str = "".join(equation_parts)[:-3]
    equation_str += f" = {next_number}"
    
    print(f"The last number is {last_number}.")
    print(f"The sum of its digits is {' + '.join(s)} = {digit_sum}.")
    print(f"The next number is calculated as: {last_number} + {digit_sum} = {next_number}.")
    print("\nFinal Equation:")
    print(equation_str)

solve()
<<<2364>>>