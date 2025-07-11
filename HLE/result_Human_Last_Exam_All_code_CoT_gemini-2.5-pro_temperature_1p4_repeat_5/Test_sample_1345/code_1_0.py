import math

def solve():
    """
    Calculates the maximal possible number of complex zeros for the given matrix determinant problem.
    """
    try:
        n_str = input("Enter the value of N: ")
        N = int(n_str)
        if N <= 0:
            print("N must be a positive integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    if N == 1:
        # For N=1, the equation is A_11 + k_1 = 0, so k_1 = -A_11.
        # Since A_11 is real, k_1 is real.
        # There are no complex zeros.
        num_zeros = 0
        calculation_str = "0"
        print("For N = 1, the equation is det(B) = A_11 + k_1 = 0.")
        print("The only solution is k_1 = -A_11, which is always real.")
    else:
        # For N > 1, the total number of solutions to the system of polynomial equations
        # is N * 2^(N-1). For a generic choice of the matrix A, all these solutions
        # correspond to complex values k_j (not purely real or imaginary).
        # The degree of the polynomial in k_1 is N * 2^(N-1), which is even for N>1.
        # A generic polynomial of even degree with real coefficients can have all its roots
        # be complex, so this maximum is attainable.
        num_zeros = N * (2**(N - 1))
        calculation_str = f"{N} * 2^({N} - 1) = {N} * {2**(N-1)} = {num_zeros}"
        print(f"For N = {N}, the maximal number of complex zeros is given by the formula N * 2^(N-1).")
        print(f"Calculation: {calculation_str}")

    print(f"\nThe maximal possible number of complex zeros is: {num_zeros}")

solve()