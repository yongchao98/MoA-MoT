import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given problem parameters.
    """
    # Given parameters
    N = 200
    k = 13

    # List to store the results
    results = []

    print(f"Calculating the number of n-simplices for N = {N}, k = {k}, and n <= 5.")
    print("The number of n-simplices corresponds to the number of sequences")
    print(f"j_0, j_1, ..., j_n such that {k} <= j_0 <= j_1 <= ... <= j_n <= {N}.")
    print("The formula is C(N - k + n + 1, n + 1).\n")

    # Loop for n from 0 to 5
    for n in range(6):
        # The number of elements in the set {k, ..., N} is N - k + 1.
        # We are choosing n+1 elements from this set with replacement.
        # The formula for combinations with repetition is C(m + r - 1, r),
        # where m = N - k + 1 and r = n + 1.
        # This simplifies to C(N - k + n + 1, n + 1).
        
        # Values for the binomial coefficient C(m_comb, r_comb)
        m_comb = N - k + n + 1
        r_comb = n + 1
        
        # Calculate the result using math.comb
        num_simplices = math.comb(m_comb, r_comb)
        results.append(num_simplices)
        
        # Print the equation with the numbers filled in
        print(f"For n = {n}:")
        print(f"Number of simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({m_comb}, {r_comb}) = {num_simplices}")
        if n < 5:
            print("-" * 30)

if __name__ == '__main__':
    solve_simplices_count()