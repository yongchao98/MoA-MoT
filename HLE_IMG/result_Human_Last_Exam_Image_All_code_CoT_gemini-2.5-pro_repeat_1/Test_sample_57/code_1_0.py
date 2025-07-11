def fibonacci(n):
    """
    Calculates the n-th Fibonacci number based on F_1=1, F_2=1.
    """
    if n <= 0:
        return 0
    elif n == 1 or n == 2:
        return 1
    
    a, b = 1, 1
    for _ in range(n - 2):
        a, b = b, a + b
    return b

def count_tiling_ways(n):
    """
    Calculates the number of ways to tile the L-shaped region of size n.
    The formula is 2 * F_n.
    """
    if n < 2:
        print(f"The shape is not well-defined for n={n} (must be n >= 2).")
        return

    fn = fibonacci(n)
    result = 2 * fn
    
    print(f"For n = {n}:")
    print(f"The number of ways is given by the formula: 2 * F_n")
    print(f"The {n}-th Fibonacci number, F_{n}, is {fn}.")
    print(f"So, the total number of ways is 2 * {fn} = {result}.")

# You can change the value of n to find the number of ways for a different size.
n = 10
count_tiling_ways(n)

# Let's also check for n=2 and n=3 to verify our initial calculations.
print("\n--- Verification ---")
count_tiling_ways(2)
print("-" * 20)
count_tiling_ways(3)