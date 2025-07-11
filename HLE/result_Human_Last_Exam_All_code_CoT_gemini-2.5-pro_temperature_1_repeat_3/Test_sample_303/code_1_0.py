def calculate_and_print_fibonacci(n):
    """
    This function calculates and prints the first n Fibonacci numbers.
    The sequence starts with 0 and 1.
    """
    if n <= 0:
        print("The number of terms must be a positive integer.")
        return

    # Initialize the first two numbers of the sequence
    a, b = 0, 1

    print(f"The first {n} Fibonacci numbers are:")

    # Loop to generate and print the numbers
    for i in range(n):
        # The equation for the sequence is F(i) = F(i-1) + F(i-2)
        # We will print each number, F(i), that is part of the final sequence.
        print(a, end=" ")
        
        # Calculate the next number and update the sequence
        next_fib = a + b
        a = b
        b = next_fib
    
    print() # Add a newline at the end for clean output

# The recipe's goal was to get the first 100 numbers.
calculate_and_print_fibonacci(100)