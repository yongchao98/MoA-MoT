def solve_sequence_count():
    """
    Calculates the number of possible sequences based on the given constraints.
    The user can modify the values of N, M, and K below.
    """
    # --- Parameters of the problem ---
    # N: The maximum value for any number in the sequence.
    # M: The maximum increase between consecutive numbers.
    # K: The length of the sequence.
    # Condition: M * (K - 1) < N must hold.
    N = 10
    M = 3
    K = 4

    print(f"Calculating the number of sequences for N={N}, M={M}, K={K}\n")

    # --- Helper function to calculate combinations C(n, k) ---
    def combinations(n, k):
        """Calculates the binomial coefficient C(n, k) = n! / (k! * (n-k)!)."""
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

    # --- Main calculation using the inclusion-exclusion formula ---
    # Formula: Sum_{j=0 to floor((N-K)/M)} [(-1)^j * C(K-1, j) * C(N - j*M, K)]
    
    total_sequences = 0
    limit = (N - K) // M

    # Build strings for printing the equation
    equation_symbolic = []
    equation_numeric = []
    equation_evaluated = []

    for j in range(limit + 1):
        comb1_n = K - 1
        comb1_k = j
        comb2_n = N - j * M
        comb2_k = K

        val_comb1 = combinations(comb1_n, comb1_k)
        val_comb2 = combinations(comb2_n, comb2_k)
        
        term_value = val_comb1 * val_comb2

        if term_value == 0:
            continue

        if j % 2 == 1:  # Odd j, subtract the term
            total_sequences -= term_value
            sign_symbol = "-"
            sign_value = f"- {term_value}"
        else:  # Even j, add the term
            total_sequences += term_value
            sign_symbol = "+"
            sign_value = f"+ {term_value}"

        # For the first term, don't print the '+' sign
        if j == 0:
            sign_symbol = ""
            sign_value = str(term_value)

        equation_symbolic.append(f"{sign_symbol} C({comb1_n}, {comb1_k}) * C({comb2_n}, {comb2_k})")
        equation_numeric.append(f"{sign_symbol} {val_comb1} * {val_comb2}")
        if j > 0:
            equation_evaluated.append(sign_value)
        else:
            equation_evaluated.append(str(term_value))


    # --- Print the results ---
    print("The number of sequences is given by the formula:")
    print(f"Result = {' '.join(equation_symbolic).lstrip('+ ')}")
    print("\nEvaluating the terms:")
    print(f"Result = {' '.join(equation_numeric).lstrip('+ ')}")
    print("\nCalculating the final sum:")
    print(f"Result = {' '.join(equation_evaluated)} = {total_sequences}")


# Execute the function
solve_sequence_count()