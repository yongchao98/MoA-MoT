def solve_sequence():
    """
    Calculates the next number in a sequence based on the rule:
    next_number = current_number + sum_of_digits(current_number).
    """
    last_number = 2352
    
    # Calculate the sum of the digits of the last number
    digits_str = str(last_number)
    digits = [int(d) for d in digits_str]
    sum_of_digits = sum(digits)
    
    # Calculate the next number
    next_number = last_number + sum_of_digits
    
    # Format the output to show the equation clearly
    digits_sum_str = " + ".join(digits_str)
    
    print(f"The last number is {last_number}.")
    print(f"The sum of its digits is {digits_sum_str} = {sum_of_digits}.")
    print("The rule to find the next number is: current number + sum of its digits.")
    print(f"So, the final equation is:")
    print(f"{last_number} + {sum_of_digits} = {next_number}")

solve_sequence()