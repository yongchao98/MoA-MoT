def solve():
    """
    This function calculates and prints the output of the Piet program.
    The Piet program computes a Fibonacci-like sequence.
    """
    # Initialize the first two numbers of the sequence.
    a, b = 1, 1

    # The program generates 6 numbers in the sequence.
    print("The Piet program performs the following calculations:")
    for _ in range(6):
        # Calculate the next number in the sequence
        result = a + b
        
        # Print the equation for the current step
        print(f"{a} + {b} = {result}")
        
        # Update the numbers for the next iteration
        a, b = b, result

solve()