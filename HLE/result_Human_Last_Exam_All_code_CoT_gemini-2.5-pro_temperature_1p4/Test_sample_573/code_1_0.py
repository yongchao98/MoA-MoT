import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given problem.
    """
    # Given parameters
    N = 200
    k = 13

    # Explain the context and derived formula
    print(f"We are finding the number of n-simplices in N_•(Z_N)_k/ for N={N} and k={k}.")
    print("This corresponds to counting the number of non-decreasing sequences of n+1 integers")
    print(f"(j_0, j_1, ..., j_n) such that {k} <= j_0 <= j_1 <= ... <= j_n <= {N}.")
    print("\nThis is a combination with repetition problem. The formula is C(m + r - 1, r),")
    print(f"where m is the size of the set of integers [{k}, {N}], so m = {N} - {k} + 1 = {N - k + 1},")
    print("and r is the length of the sequence, so r = n + 1.")
    print("The final formula is C(N - k + n + 1, n + 1).\n")
    print("-" * 50)

    # Loop for n from 0 to 5
    for n in range(6):
        # Parameters for the binomial coefficient C(p, q)
        p = N - k + n + 1
        q = n + 1

        # Calculate the number of simplices using math.comb
        try:
            num_simplices = math.comb(p, q)
            print(f"For n = {n}:")
            # Display the full calculation expression as requested
            print(f"Number of {n}-simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({p}, {q}) = {num_simplices}")
        except ValueError as e:
            print(f"Could not compute for n = {n}. Error: {e}")
        
        print("-" * 50)

# Run the solver
solve_simplices_count()