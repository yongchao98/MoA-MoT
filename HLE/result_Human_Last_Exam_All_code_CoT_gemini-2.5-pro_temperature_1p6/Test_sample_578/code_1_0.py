import math

def catalan(n):
    """Computes the n-th Catalan number."""
    if n < 0:
        return 0
    # Using the formula C_n = comb(2n, n) // (n + 1)
    return math.comb(2 * n, n) // (n + 1)

def solve():
    """
    Calculates the product of d_n for n from 1 to 8.
    d_n is the dimension of the B_n-fixed subspace of Kh(T(n,n); Q),
    which is hypothesized to be 2 * C_n.
    """
    product = 1
    d_values = []
    
    for n in range(1, 9):
        c_n = catalan(n)
        d_n = 2 * c_n
        d_values.append(d_n)
        product *= d_n
        
    # Formatting the output string with the full equation.
    d_str = " * ".join(map(str, d_values))
    print(f"The values for d_n from n=1 to 8 are: {d_values}")
    print("The product is d_1 * d_2 * d_3 * d_4 * d_5 * d_6 * d_7 * d_8")
    print(f"= {d_str}")
    print(f"= {product}")

solve()