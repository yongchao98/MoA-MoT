import sys

# The Chef script is an allegory for a recursive function. The "Caramel Sauce"
# recipe calls itself, which is recursion. However, it lacks a "base case"—a
# condition to stop recursing—so it would call itself infinitely. This is a fatal
# error. The correct fix from the choices is F, which removes the recursive call,
# thus stopping the program from crashing.
#
# The original intent was to calculate Fibonacci numbers. The phrase "an awful lot
# of caramel sauce" alludes to the inefficiency of naive recursive Fibonacci, which
# re-computes the same values many times.
#
# This Python code fulfills the original intent of the script correctly and efficiently
# by using memoization (a cache) to avoid re-computing values.

# The original recipe mentioned calculating up to 100, which might exceed the
# default recursion limit in a simple recursive solution.
sys.setrecursionlimit(2000)

# Memoization cache to store already computed Fibonacci numbers
memo = {}

def fibonacci(n):
    """
    Calculates the nth Fibonacci number using recursion with memoization
    to make it efficient.
    """
    if n in memo:
        return memo[n]
    if n < 0:
        raise ValueError("Input must be a non-negative integer.")
    # Base cases for the recursion
    if n <= 1:
        return n

    # The recursive step: F(n) = F(n-1) + F(n-2)
    result = fibonacci(n - 1) + fibonacci(n - 2)

    # Cache the result before returning
    memo[n] = result
    return result

# Let's calculate the 10th Fibonacci number as an example.
n = 10

# To fulfill the requirement "output each number in the final equation!",
# we will show the final step of the calculation for F(10).
if n > 1:
    # We calculate the final number, and the two numbers that sum to it.
    # Our efficient function uses the cache, so these calls are very fast.
    final_result = fibonacci(n)
    n_minus_1 = memo.get(n - 1, fibonacci(n - 1)) # Get from cache if available
    n_minus_2 = memo.get(n - 2, fibonacci(n - 2)) # Get from cache if available

    print(f"The equation for Fibonacci({n}) is: F({n}) = F({n-1}) + F({n-2})")
    print("The final calculation with each number is:")
    print(f"{final_result} = {n_minus_1} + {n_minus_2}")
else:
    # Handle the base cases 0 and 1
    final_result = fibonacci(n)
    print(f"Fibonacci({n}) = {final_result}")
