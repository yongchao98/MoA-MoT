def find_next_harshad_number():
    """
    This function finds the next number in a sequence of Harshad numbers.
    A Harshad number is an integer that is divisible by the sum of its digits.
    """
    
    # The last number in the provided sequence
    last_number = 2352
    
    # Start checking from the next integer
    num = last_number + 1
    
    next_harshad = 0
    
    while True:
        # Calculate the sum of the digits of the current number
        s = sum(int(digit) for digit in str(num))
        
        # Check if the number is divisible by the sum of its digits
        if s != 0 and num % s == 0:
            next_harshad = num
            break
        
        num += 1

    # Explain the pattern and the solution
    print("The pattern in the sequence is that each number is a Harshad number.")
    print("A Harshad number is an integer that is divisible by the sum of its digits.")
    print(f"The last number given is {last_number}.")
    print(f"To find the next number, we search for the smallest integer greater than {last_number} that is a Harshad number.")
    print(f"The next Harshad number found is {next_harshad}.")
    
    print("\nHere is the verification for the next number:")
    
    digits_str = list(str(next_harshad))
    sum_of_digits = sum(int(d) for d in digits_str)
    
    print(f"The number is {next_harshad}.")
    print(f"The sum of its digits is: {' + '.join(digits_str)} = {sum_of_digits}")
    
    # Display the final equation as requested
    print("\nThe final equation showing it's a Harshad number is:")
    print(f"{next_harshad} % ({' + '.join(digits_str)}) = {next_harshad % sum_of_digits}")
    
    print(f"\nSince the remainder is 0, {next_harshad} is the next number in the sequence.")

# Execute the function to find and print the result
find_next_harshad_number()