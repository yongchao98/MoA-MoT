import math

def catalan(n):
    """Computes the n-th Catalan number."""
    if n < 0:
        return 0
    # C_n = (1/(n+1)) * (2n choose n)
    return math.comb(2 * n, n) // (n + 1)

def get_d(n):
    """
    Computes the dimension d_n of the B_n-invariant subspace of Kh(T(n,n); Q).
    """
    if n == 1:
        # For n=1, T(1,1) is the unknot. B_1 is trivial, so the whole
        # 2-dimensional Khovanov homology is the invariant subspace.
        return 2
    else:
        # For n>=2, d_n is given by the n-th Catalan number.
        return catalan(n)

def solve():
    """
    Calculates the product of d_n for n from 1 to 8.
    """
    d_values = []
    for n in range(1, 9):
        d_n = get_d(n)
        d_values.append(d_n)
        print(f"d_{n} = {d_n}")

    product = 1
    for val in d_values:
        product *= val
    
    equation_str = " * ".join(map(str, d_values))
    print(f"\nThe product is {equation_str} = {product}")

solve()