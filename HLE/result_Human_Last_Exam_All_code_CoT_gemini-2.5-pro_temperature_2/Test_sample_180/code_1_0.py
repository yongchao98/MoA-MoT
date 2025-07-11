import fractions

def find_lambda_candidate():
    """
    This function searches for the value of lambda by calculating the ratio
    |Q(A)| / |A|^4 for sets A_n = {2^0, 2^1, ..., 2^(n-1)}. This family of
    sets is known to generate a large number of distinct differences, making
    them good candidates for maximizing the ratio.
    """
    max_ratio = 0
    lambda_candidate = 0
    final_eq = ""

    print("Investigating the ratio |Q(A)| / |A|^4 for sets A_n = {2^0, 2^1, ..., 2^(n-1)}:")
    print("-" * 70)

    for n in range(2, 6):
        # Create the set A = {1, 2, 4, ..., 2^(n-1)}
        A = {2**i for i in range(n)}
        
        # Calculate the set of differences D(A) = {a-b | a,b in A}
        diffs = set()
        for a in A:
            for b in A:
                diffs.add(a - b)
        
        # The set of non-zero differences for the denominator
        non_zero_diffs = {d for d in diffs if d != 0}

        # Calculate the quotient set Q(A)
        # Using Fraction for exact rational arithmetic
        quotients = {fractions.Fraction(0, 1)}
        for a_minus_b in diffs:
            for c_minus_d in non_zero_diffs:
                 quotients.add(fractions.Fraction(a_minus_b, c_minus_d))
        
        q_size = len(quotients)
        n_4 = n**4
        ratio = q_size / n_4
        
        print(f"For n={n}, A={sorted(list(A))}:")
        print(f"|Q(A)| = {q_size}")
        print(f"|A|^4 = {n}^4 = {n_4}")
        print(f"Ratio: {q_size}/{n_4} \u2248 {ratio:.6f}")
        print("-" * 25)

        if ratio > max_ratio:
            max_ratio = ratio
            lambda_candidate = fractions.Fraction(q_size, n_4)
            final_eq = f"{q_size} / {n}^4 = {q_size} / {n_4}"
    
    print("\nThe exploration suggests the maximum ratio is achieved at n=2.")
    print(f"The best candidate for lambda is {lambda_candidate}, which is {float(lambda_candidate)}.")
    print(f"The equation giving this value is |Q(A)| / |A|^4 = {final_eq}.")

find_lambda_candidate()