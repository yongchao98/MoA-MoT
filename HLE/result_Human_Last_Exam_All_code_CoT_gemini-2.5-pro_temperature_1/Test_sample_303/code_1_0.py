def print_fibonacci_sequence(count):
    """
    Calculates and prints the first 'count' numbers in the Fibonacci sequence.
    This demonstrates a correct, non-recursive way to generate the sequence,
    contrasting with the flawed recursive logic in the Chef script.
    """
    print(f"The first {count} Fibonacci numbers are:")
    
    # Initialize the first two numbers
    a, b = 0, 1
    
    # Handle the first two numbers being printed
    if count >= 1:
        print(a, end=" ")
    if count >= 2:
        print(b, end=" ")
    
    # Generate and print the rest of the sequence
    for i in range(2, count):
        next_fib = a + b
        print(next_fib, end=" ")
        # Update the sequence variables
        a, b = b, next_fib
        
    print() # for a final newline

# The original script mentions 100 numbers, but we will demonstrate with the first 10.
print_fibonacci_sequence(10)