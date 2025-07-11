def solve_fibonacci():
    """
    This function calculates and prints the first 100 Fibonacci numbers.

    The "Chef script" puzzle contains a logical error of infinite recursion. The
    "Caramel Sauce" recipe calls itself without a base case to stop, leading to
    an endless loop (metaphorically, "an awful lot of caramel sauce"). The correct
    fix from the options is F, which removes the self-referential call, thus
    stopping the infinite loop.

    This script correctly implements the intended logic of the recipe. It
    calculates the sequence iteratively, which is efficient, and prints each
    of the first 100 numbers (from F(0) to F(99)) in an equation format.
    """
    
    # Initialize the first two numbers of the sequence
    a, b = 0, 1
    
    # The first 100 Fibonacci numbers corresponds to n from 0 to 99
    limit = 100
    
    print(f"Printing the first {limit} Fibonacci numbers:")

    for n in range(limit):
        if n == 0:
            # The first number in the sequence
            result = a
            # The "equation" for the base case F(0)
            print(f"F(0) = {result}")
        elif n == 1:
            # The second number in the sequence
            result = b
            # The "equation" for the base case F(1)
            print(f"F(1) = {result}")
        else:
            # Calculate the next number in the sequence
            next_fib = a + b
            # The "equation" shows the resulting number
            print(f"F({n}) = {next_fib}")
            
            # Update the last two numbers for the next iteration
            a = b
            b = next_fib

solve_fibonacci()