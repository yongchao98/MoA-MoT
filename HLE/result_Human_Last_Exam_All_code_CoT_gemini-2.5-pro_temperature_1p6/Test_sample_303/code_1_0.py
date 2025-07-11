def print_fibonacci_sequence(count):
    """
    This function calculates and prints a specified number of Fibonacci numbers.

    The Fibonacci sequence starts with 0 and 1, and each subsequent number is the
    sum of the two preceding ones (e.g., 0, 1, 1, 2, 3, 5, ...).

    Args:
        count (int): The number of Fibonacci numbers to print.
    """
    a, b = 0, 1  # Initialize the first two numbers in the sequence
    
    # Check if the requested count is valid
    if count <= 0:
        print("Please enter a positive integer for the count.")
        return
        
    print(f"The first {count} Fibonacci numbers are:")
    # Loop 'count' times to generate and print each number
    for i in range(count):
        print(f"F({i}): {a}")
        # Update the values for the next iteration in the sequence
        a, b = b, a + b

# The original script intended to print the first 100 Fibonacci numbers.
print_fibonacci_sequence(100)