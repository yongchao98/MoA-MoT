import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    # math.comb is available in Python 3.8+ and is recommended
    try:
        return math.comb(n, k)
    except AttributeError:
        # Fallback for older Python versions
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

def solve_sequence_count():
    """
    Prompts the user for N, M, and K, and calculates the number of
    possible sequences based on the derived formula.
    """
    try:
        # Read inputs from the user
        N = int(input("Enter the value for N (maximum value in the sequence): "))
        M = int(input("Enter the value for M (maximum increase between consecutive numbers): "))
        K = int(input("Enter the value for K (length of the sequence): "))

        # Validate the inputs based on problem constraints
        if not (N > 0 and M > 0 and K > 0):
            print("\nError: N, M, and K must be positive integers.")
            return

        # For a sequence of length K>1, the condition M(K-1) < N is specified.
        if K > 1 and not (M * (K - 1) < N):
            print(f"\nError: The condition M*(K-1) < N is not met.")
            print(f"({M}*({K}-1) = {M*(K-1)} is not strictly less than {N})")
            return
        
        # The formula is Sum_{j=0}^{K-1} (-1)^j * C(K-1, j) * C(N - M*j, K)
        total_count = 0
        
        # Build string representations for the calculation steps
        equation_symbolic = "Count = "
        equation_values = "      = "
        equation_terms = "      = "
        
        for j in range(K):
            # Calculate components of the term
            comb1 = combinations(K - 1, j)
            n_val = N - M * j
            comb2 = combinations(n_val, K)
            term_value = comb1 * comb2
            
            # Append to the symbolic formula string
            if j > 0:
                equation_symbolic += " - " if j % 2 == 1 else " + "
            equation_symbolic += f"C({K-1}, {j}) * C({n_val}, {K})"
            
            # Append to the formula with combination values
            if j > 0:
                equation_values += " - " if j % 2 == 1 else " + "
            equation_values += f"{comb1}*{comb2}"
            
            # Append to the string with calculated term values
            if j > 0:
                equation_terms += " - " if j % 2 == 1 else " + "
            equation_terms += str(term_value)
            
            # Add or subtract the term from the total
            if j % 2 == 1:
                total_count -= term_value
            else:
                total_count += term_value

            # Optimization: if the upper value of a combination is less than the lower,
            # all subsequent terms will be zero.
            if n_val < K:
                break
        
        final_calculation = f"{equation_terms} = {total_count}"
        
        print("\nThe number of sequences is calculated using the formula derived from the Principle of Inclusion-Exclusion:")
        print(equation_symbolic)
        print(equation_values)
        print(final_calculation)
        print(f"\nResult: The total number of possible sequences is {total_count}.")

    except ValueError:
        print("\nError: Invalid input. Please enter valid integers for N, M, and K.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == '__main__':
    solve_sequence_count()