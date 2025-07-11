def solve_fibonacci_puzzle():
    """
    This function implements the logic described in the "recipe" to generate
    the first 100 Fibonacci numbers. A naive recursive approach, like the one
    implied by the recipe, is very inefficient (the "awful lot of caramel sauce").
    This iterative approach is much more efficient.

    The code will print the equation for each number in the sequence,
    showing how it is the sum of the previous two.
    """
    limit = 100
    
    # Initialize the first two numbers of the sequence
    a, b = 0, 1
    
    # Store the sequence of numbers
    fib_sequence = []
    
    if limit >= 1:
        fib_sequence.append(a)
    if limit >= 2:
        fib_sequence.append(b)
    
    # Generate the rest of the sequence up to the limit
    for _ in range(2, limit):
        next_fib = a + b
        fib_sequence.append(next_fib)
        a = b
        b = next_fib

    print(f"The first {limit} Fibonacci numbers and their relationships:")
    
    # Print the base cases
    if limit >= 1:
        print(f"F(0) = {fib_sequence[0]}")
    if limit >= 2:
        print(f"F(1) = {fib_sequence[1]}")
    
    # Print the equation for the remaining numbers
    # The prompt asks to "output each number in the final equation"
    for i in range(2, limit):
        current_num = fib_sequence[i]
        prev_num1 = fib_sequence[i-1]
        prev_num2 = fib_sequence[i-2]
        
        # F(n) = F(n-1) + F(n-2)
        print(f"{current_num} = {prev_num1} + {prev_num2}")

solve_fibonacci_puzzle()