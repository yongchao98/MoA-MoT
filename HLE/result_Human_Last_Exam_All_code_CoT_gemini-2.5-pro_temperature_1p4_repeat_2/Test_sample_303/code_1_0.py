def fibonacci(n):
    """
    Calculates the nth Fibonacci number using an iterative approach,
    which is efficient for large numbers.
    """
    if n < 0:
        raise ValueError("Input cannot be negative.")
    
    # Base cases
    if n == 0:
        return 0
    if n == 1:
        return 1

    # Iterative calculation
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

# The recipe aims to compute the 100th Fibonacci number.
# The recursive definition is F(100) = F(99) + F(98).
n = 100
fib_100 = fibonacci(n)
fib_99 = fibonacci(n - 1)
fib_98 = fibonacci(n - 2)

# Print each number in the final equation
print(f"{fib_100} = {fib_99} + {fib_98}")