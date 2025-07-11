def solve_task():
    """
    Calculates the maximal possible number of complex zeros for the given matrix system
    of dimension N.
    """
    try:
        # We demonstrate the calculation for a few values of N, as an example.
        # A specific N value is not given in the problem, so we present the general solution.
        n_values = [1, 2, 3, 4]
        print("Calculating the maximal number of complex zeros for N = 1, 2, 3, 4.\n")
        
        for N in n_values:
            print(f"--- For N = {N} ---")
            if N == 1:
                max_zeros = 0
                degree = N * (2**(N - 1))
                print(f"The degree of the characteristic polynomial in k1 is N * 2^(N-1) = {N} * 2^({N}-1) = {degree}.")
                print("A polynomial with real coefficients and an odd degree must have at least one real root.")
                print("This means k1 must be real, violating the complex zero condition.")
                print(f"Maximal number of complex zeros: {max_zeros}")
                print(f"Final Equation: {max_zeros}")
            else:
                # For N > 1, the degree is N * 2^(N-1) which is always even.
                # A real polynomial of even degree can have all its roots be complex.
                max_zeros = N * (2**(N - 1))
                print(f"The degree of the characteristic polynomial in k1 is N * 2^(N-1).")
                print(f"The degree is even, so all roots can be complex.")
                print("The maximal number of complex zeros is equal to this degree.")
                print(f"Maximal number of complex zeros: {N} * (2**({N}-1)) = {max_zeros}")
            
            print("-" * (17 + len(str(N))))
            print() # Add a newline for better readability
            
    except Exception as e:
        print(f"An error occurred: {e}")

solve_task()