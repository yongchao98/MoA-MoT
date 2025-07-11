def solve():
    """
    Calculates the maximal possible number of complex zeros for the determinant
    of the matrix B(k) for a given size N based on the derived formula.
    """
    # Set the size of the matrix N here.
    N = 4

    if not isinstance(N, int) or N < 0:
        print("N must be a non-negative integer.")
        return
    
    print(f"The maximal possible number of complex zeros for an N x N matrix is given by the formula d(N) = N * 2^(N-1).")
    print(f"This formula gives the degree of the polynomial in k_1 whose roots are the solutions.")
    print(f"\nCalculating for N = {N}:")

    if N == 0:
        result = 0
        power_of_2 = "N/A"
    else:
        # The formula for the maximal number of zeros is N * 2^(N-1).
        power_of_2 = 2**(N - 1)
        result = N * power_of_2

    if N > 0:
        print(f"d({N}) = {N} * 2^({N} - 1) = {N} * {power_of_2} = {result}")
    else:
        print(f"d(0) = {result}")

solve()