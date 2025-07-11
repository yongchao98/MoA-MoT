import math

def solve_problem():
    """
    Calculates the product of d_n for n from 1 to 8.
    d_n is the dimension of the B_n-fixed subspace of the Khovanov homology
    of the torus link T(n,n). It is known that d_n = n + 1.
    The problem reduces to calculating the product of (n+1) for n=1 to 8,
    which is 2 * 3 * ... * 9 = 9!.
    """

    limit = 8
    
    # Calculate each d_n and store it in a list.
    d_values = []
    for n in range(1, limit + 1):
        d_n = n + 1
        d_values.append(d_n)

    # Calculate the product of all d_n values.
    # This is equivalent to (limit + 1)! or math.factorial(9).
    product_result = math.factorial(limit + 1)

    # Create a string representation of the equation.
    equation_terms_str = " * ".join(map(str, d_values))
    
    # Print the final equation with all the numbers.
    print(f"The product is d_1 * d_2 * ... * d_8")
    print(f"Based on the formula d_n = n + 1, we need to compute:")
    print(f"{equation_terms_str} = {product_result}")

solve_problem()