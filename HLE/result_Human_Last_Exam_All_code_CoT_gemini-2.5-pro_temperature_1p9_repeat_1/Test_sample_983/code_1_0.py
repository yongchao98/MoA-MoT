def find_next_number_in_sequence():
    """
    Finds the next number in the sequence starting from after 2352.
    The rule is that the number must be divisible by 6 and also by the sum of its digits.
    """
    
    # Start checking from the number right after the last one in the sequence.
    num = 2352

    while True:
        num += 1
        
        # Rule 1: Must be divisible by 6.
        if num % 6 == 0:
            
            # Rule 2: Must be a Harshad number (divisible by the sum of its digits).
            sum_of_digits = sum(int(digit) for digit in str(num))
            
            if sum_of_digits > 0 and num % sum_of_digits == 0:
                # This is the next number in the sequence.
                print(f"The next number is {num}.")
                print(f"It satisfies the rule based on the following checks:")
                print(f"1. Check for divisibility by 6:")
                print(f"{num} % 6 = {num % 6}")
                print(f"2. Check for being a Harshad number (divisibility by sum of its digits):")
                print(f"Sum of digits of {num} is {sum_of_digits}.")
                print(f"{num} % {sum_of_digits} = {num % sum_of_digits}")
                
                return num

# Run the function to get the result.
# The value is not printed directly but returned and then printed by the final line.
result = find_next_number_in_sequence()