import math

def calculate_product_of_d_n():
    """
    Calculates the product of d_n for n=1 to 8, where d_n is the
    dimension of the B_n-invariant subspace of Kh(T(n,n); Q).

    The dimension d_n is given by the (n+2)-th Fibonacci number, F_{n+2}.
    """

    # Step 1: Define a function to compute Fibonacci numbers.
    # We use the standard definition F_1=1, F_2=1.
    def fibonacci(k):
        if k <= 0:
            return 0
        elif k == 1:
            return 1
        
        a, b = 1, 1
        for _ in range(k - 2):
            a, b = b, a + b
        return b

    # Step 2: Generate the list of d_n values.
    # d_n = F_{n+2} for n=1, ..., 8.
    d_values = []
    for n in range(1, 9):
        # We need F_{n+2}.
        fib_index = n + 2
        d_n = fibonacci(fib_index)
        d_values.append(d_n)

    # Step 3: Calculate the product.
    product = math.prod(d_values)

    # Step 4: Print the output in the desired format.
    equation_str = " * ".join(map(str, d_values))
    print(f"Based on the formula d_n = F_(n+2), the values for n=1 to 8 are:")
    for i, val in enumerate(d_values, 1):
        print(f"d_{i} = F_{i+2} = {val}")

    print("\nThe product is calculated as:")
    print(f"{equation_str} = {product}")

calculate_product_of_d_n()