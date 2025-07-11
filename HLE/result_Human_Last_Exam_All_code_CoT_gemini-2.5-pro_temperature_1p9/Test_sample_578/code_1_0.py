def fibonacci(k):
    """
    Calculates the k-th Fibonacci number, with the standard
    convention F_1 = 1, F_2 = 1.
    """
    if k <= 0:
        raise ValueError("Input must be a positive integer.")
    
    # F_1 = 1, F_2 = 1
    if k == 1:
        return 1
    if k == 2:
        return 1

    a, b = 1, 1
    # We need k-2 iterations to get from F_2 to F_k
    for _ in range(k - 2):
        a, b = b, a + b
    return b

def solve():
    """
    Calculates the product of d_n for n from 1 to 8.
    """
    d_values = []
    for n in range(1, 9):
        # The dimension d_n is given by the (n+2)-th Fibonacci number.
        fib_index = n + 2
        d_n = fibonacci(fib_index)
        d_values.append(d_n)

    product = 1
    for value in d_values:
        product *= value
    
    # Construct the output string showing the full equation.
    equation_str = " * ".join(map(str, d_values))
    
    print("The values for d_n, the dimension of the B_n-invariant subspace of Kh(T(n,n); Q), are given by Fibonacci numbers F_{n+2}.")
    print("For n = 1 to 8, the values are:")
    for i in range(len(d_values)):
        print(f"d_{i+1} = F_{i+3} = {d_values[i]}")
    
    print("\nThe product is:")
    print(f"{equation_str} = {product}")

solve()