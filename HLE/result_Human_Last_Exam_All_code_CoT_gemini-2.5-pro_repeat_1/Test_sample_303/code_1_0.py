def print_fibonacci_sequence(count):
    """
    Calculates and prints a specified number of Fibonacci numbers.
    The original Chef script was intended to do this for the first 100 numbers.
    """
    # Initialize the first two numbers of the sequence
    num1, num2 = 0, 1
    
    if count <= 0:
        print("Please provide a positive number.")
        return
        
    print(f"The first {count} Fibonacci numbers:")
    
    # Loop to generate and print the sequence
    for i in range(count):
        # Print the current number in the sequence
        print(num1)
        
        # The equation for the next number is: next_num = num1 + num2
        # Here, we are outputting `num1`, which is one of the numbers
        # in the equation for the next step.
        next_num = num1 + num2
        
        # Update the numbers for the next iteration
        num1 = num2
        num2 = next_num

# Execute the function to print the first 100 Fibonacci numbers
print_fibonacci_sequence(100)