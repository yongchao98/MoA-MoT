import math

def calculate_max_complex_zeros(N):
    """
    Calculates the maximal possible number of complex zeros for the given NxN matrix B(k).

    Args:
        N (int): The dimension of the matrix.

    Returns:
        int: The maximal number of complex zeros.
    """
    if not isinstance(N, int) or N <= 0:
        raise ValueError("N must be a positive integer.")

    if N == 1:
        # For N=1, the equation is linear with a real solution, so 0 complex zeros.
        return 0
    else:
        # For N > 1, the number of zeros is given by the degree of the characteristic polynomial,
        # which is N * 2^(N-1). For a generic choice of matrix A, all these zeros can be complex.
        return N * (2**(N - 1))

def main():
    """
    Main function to demonstrate the calculation for various N.
    """
    print("This program calculates the maximal possible number of complex zeros for the matrix B(k) of size N x N.")
    print("The formula is 0 for N=1, and N * 2^(N-1) for N > 1.\n")
    
    test_cases = [1, 2, 3, 4, 5, 10]
    
    for n in test_cases:
        result = calculate_max_complex_zeros(n)
        print(f"For N = {n}:")
        if n == 1:
            print(f"Maximal number of complex zeros = {result}")
        else:
            # Outputting each number in the final equation
            power = n - 1
            base = 2
            multiplier = n
            print(f"Maximal number of complex zeros = {multiplier} * {base}^{power} = {result}")
        print("-" * 20)

if __name__ == "__main__":
    main()
