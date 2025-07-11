def fib(n):
    """
    Calculates the nth Fibonacci number using an iterative approach.
    """
    if n < 0:
        raise ValueError("Input must be a non-negative integer.")
    if n == 0:
        return 0
    
    # Start with the first two numbers in the sequence
    a, b = 0, 1
    # Loop n-1 times to get to the nth number
    for _ in range(n - 1):
        a, b = b, a + b
    return b

# The Chef script's main ingredient "100 g flour" implies it's trying to
# calculate the 100th Fibonacci number.
n = 100

# Calculate the required Fibonacci numbers for the equation.
fib_n = fib(n)
fib_n_minus_1 = fib(n - 1)
fib_n_minus_2 = fib(n - 2)

# Print the final equation as requested.
print(f"The error in the Chef script is an infinite recursion.")
print(f"The fix is to remove the first recursive call, 'Serve with caramel sauce.'")
print(f"\nA corrected program would calculate fib({n}) using the equation fib(n) = fib(n-1) + fib(n-2).")
print(f"\nThe final equation for n={n} is:")
print(f"{fib_n_minus_1} + {fib_n_minus_2} = {fib_n}")