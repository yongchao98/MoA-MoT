def fibonacci(n):
    """
    Calculates the n-th Fibonacci number based on F_1=1, F_2=1.
    """
    if n <= 0:
        return 0
    if n == 1 or n == 2:
        return 1
    
    a, b = 1, 1
    for _ in range(3, n + 1):
        a, b = b, a + b
    return b

def solve_tiling(n):
    """
    Calculates the number of ways to tile the L-shaped region.
    """
    if n < 2:
        print(f"The shape is not well-defined for n={n} (must be n>=2).")
        return

    # The number of ways is given by the formula 2 * F_n
    fn = fibonacci(n)
    total_ways = 2 * fn
    
    print(f"The number of ways to tile the shape for a given n is 2 * F_n.")
    print(f"The problem defines the Fibonacci sequence as F_1=1, F_2=1, F_3=2, ...")
    print(f"\nAs an example, let's calculate the number of ways for n = {n}:")
    print(f"First, we find F_{n} which is the {n}-th Fibonacci number.")
    print(f"F_{n} = {fn}")
    print(f"Then, the total number of ways is 2 * F_{n}.")
    print(f"So, the calculation is: 2 * {fn} = {total_ways}")
    print(f"\nThus, for n={n}, there are {total_ways} ways to tile the shape.")

# We can demonstrate the solution for a sample value, e.g., n=10.
# The user can change this value to calculate for other n.
n_example = 10
solve_tiling(n_example)