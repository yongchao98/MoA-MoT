import math

def solve_problem():
    """
    Calculates the product of the dimensions of B_n-invariant subspaces of Khovanov homology of T(n,n).
    """

    def catalan(m):
        """
        Calculates the m-th Catalan number, C_m.
        """
        if m < 0:
            return 0
        # C_m = (1/(m+1)) * (2m choose m)
        return math.comb(2 * m, m) // (m + 1)

    def get_d(n):
        """
        Calculates d_n, the dimension of the B_n-invariant subspace of Kh(T(n,n); Q).
        The formula is d_n = 2 * C_{n-1}.
        """
        if n <= 0:
            raise ValueError("n must be a positive integer.")
        return 2 * catalan(n - 1)

    print("The problem is to compute the product of d_n for n from 1 to 8.")
    print("The dimension d_n is given by the formula d_n = 2 * C_{n-1}, where C_m is the m-th Catalan number.")
    print("-" * 30)
    print("Calculating the values of d_n:")

    d_values = []
    product = 1

    for n in range(1, 9):
        val = get_d(n)
        d_values.append(val)
        product *= val
        print(f"For n={n}, d_{n} = 2 * C_{n-1} = 2 * C_{n-1} = {val}")

    print("-" * 30)
    print("Now, we compute the product of these values.")
    
    equation_str = " * ".join(map(str, d_values))
    
    print("The final equation is:")
    print(f"{equation_str} = {product}")

solve_problem()