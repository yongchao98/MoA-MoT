def print_fibonacci_sequence(count):
    """
    Calculates and prints the first 'count' numbers of the Fibonacci sequence.
    The sequence is the "final equation" requested, and each number is outputted.
    """
    a, b = 0, 1
    print(f"The first {count} Fibonacci numbers are:")
    if count <= 0:
        print("Please provide a positive count.")
        return
        
    for i in range(count):
        print(a, end=" ")
        # This is the core logic of the Fibonacci sequence: F(n) = F(n-1) + F(n-2)
        # We calculate the next number and then update our two pointers.
        a, b = b, a + b
    
    # Print a newline at the end for better formatting
    print()

# The original script intended to print the first 100 Fibonacci numbers.
print_fibonacci_sequence(100)