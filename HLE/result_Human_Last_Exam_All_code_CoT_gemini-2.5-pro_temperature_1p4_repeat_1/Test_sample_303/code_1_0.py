def solve_fibonacci():
    """
    This function correctly calculates and prints the first 100 Fibonacci numbers.
    Unlike the Chef recipe, this implementation is iterative and efficient,
    and it correctly terminates.
    """
    n = 100
    # Initialize the first two numbers of the sequence
    a, b = 0, 1
    
    if n <= 0:
        print("Number of terms must be a positive integer.")
        return

    print("The first 100 Fibonacci numbers are:")
    
    # We will print each number in the resulting sequence.
    output = []
    for _ in range(n):
        output.append(str(a))
        # This is the core equation: F(n) = F(n-1) + F(n-2)
        # We calculate the next number and update the previous two.
        next_val = a + b
        a = b
        b = next_val
        
    print(", ".join(output))

solve_fibonacci()