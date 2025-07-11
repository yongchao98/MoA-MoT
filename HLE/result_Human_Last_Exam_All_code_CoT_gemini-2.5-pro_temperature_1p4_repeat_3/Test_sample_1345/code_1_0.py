def calculate_max_complex_zeros():
    """
    Calculates the maximal possible number of complex zeros for the given matrix determinant problem.
    The user is prompted to enter the size of the matrix, N.
    """
    try:
        n_str = input("Enter the matrix size N (a positive integer): ")
        N = int(n_str)
        if N < 1:
            print("Error: N must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    print("-" * 20)
    if N == 1:
        # For N=1, the degree D = 1 * 2^(1-1) = 1.
        # A real polynomial of odd degree must have at least one real root.
        # The maximum number of complex roots for k_1^2 is D-1=0.
        degree_D = 1
        max_zeros = 0
        print(f"For N = {N}:")
        print(f"The degree of the polynomial in x=k_1^2 is D = {N} * 2^({N}-1) = {degree_D}.")
        print("Since D is odd, the polynomial must have at least one real root.")
        print("The maximum number of complex roots for x is D - 1 = 0.")
        print(f"Thus, the maximal number of complex zeros is {max_zeros}.")

    else: # N >= 2
        # For N>=2, the degree D = N * 2^(N-1) is always even.
        # A real polynomial of even degree can have all its roots complex.
        # Each complex root x for k_1^2 gives two complex zeros for k_1.
        # The maximal number is 2 * D = 2 * N * 2^(N-1) = N * 2^N.
        power_D = N - 1
        term_2_D = 2**power_D
        degree_D = N * term_2_D
        
        power_zeros = N
        term_2_zeros = 2**power_zeros
        max_zeros = N * term_2_zeros
        
        print(f"For N = {N}:")
        print(f"The degree of the polynomial in x=k_1^2 is D = {N} * 2^({N}-1) = {N} * {term_2_D} = {degree_D}.")
        print("Since D is even, it's possible for all D roots to be complex.")
        print("Each complex root for x gives 2 complex zeros for k_1.")
        print("The maximal number of complex zeros is 2 * D = N * 2^N.")
        print(f"Calculation: {N} * 2^{power_zeros} = {N} * {term_2_zeros} = {max_zeros}.")

if __name__ == '__main__':
    calculate_max_complex_zeros()