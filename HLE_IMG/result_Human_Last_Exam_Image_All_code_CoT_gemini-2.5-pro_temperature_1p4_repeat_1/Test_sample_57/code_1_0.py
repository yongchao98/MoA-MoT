def fibonacci(k):
    """Calculates the k-th Fibonacci number (F_1=1, F_2=1)."""
    if k <= 0:
        return 0
    if k <= 2:
        return 1
    a, b = 1, 1
    for _ in range(3, k + 1):
        a, b = b, a + b
    return b

def calculate_tiling_ways(n):
    """
    Calculates the number of ways to tile the L-shaped region of size n.
    """
    if n < 2:
        print("The shape is not well-defined for n < 2.")
        return

    # Calculate the required Fibonacci numbers
    fn = fibonacci(n)
    fn_minus_1 = fibonacci(n - 1)

    # Calculate the total number of ways using the formula g(n) = Fn^2 + F_{n-1}^2
    total_ways = fn**2 + fn_minus_1**2

    # Print the result in the format of the equation
    print(f"For n = {n}:")
    print(f"The number of ways is given by the formula F_n^2 + F_{n-1}^2.")
    print(f"Substituting the values, we get:")
    print(f"({fn})^2 + ({fn_minus_1})^2 = {fn**2} + {fn_minus_1**2} = {total_ways}")

# You can change this value to calculate for a different n
n_example = 10
calculate_tiling_ways(n_example)

# Let's also verify using the other formula F_{2n-1}
ways_alternate = fibonacci(2 * n_example - 1)
print(f"\nVerification using the formula F_(2n-1):")
print(f"F_({2 * n_example - 1}) = F_{2*n_example-1} = {ways_alternate}")
if (fibonacci(n_example)**2 + fibonacci(n_example - 1)**2) == ways_alternate:
    print("The results match.")
else:
    print("The results do not match.")
