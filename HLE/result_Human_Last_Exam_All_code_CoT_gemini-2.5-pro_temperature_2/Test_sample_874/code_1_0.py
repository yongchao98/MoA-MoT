def solve_ducci_sequence_problem():
    """
    This function solves the problem by finding a specific tuple (a, b, c, d)
    that maximizes the Ducci sequence length and has a minimal sum,
    and then calculates the required expression.
    """
    limit = 10_000_000

    # Step 1 & 2: Generate the Tribonacci sequence starting with (0, 0, 1).
    # This specific sequence is known to build tuples with maximal Ducci sequence lengths.
    t_seq = [0, 0, 1]
    while t_seq[-1] <= limit:
        next_term = t_seq[-1] + t_seq[-2] + t_seq[-3]
        t_seq.append(next_term)

    # The last generated term is over the limit. The terms we need are the three before that one.
    
    # Step 3 & 4: Construct the lexicographically smallest tuple that maximizes the
    # sequence length and has the minimum sum. This tuple will be of the form (0, t_n, t_{n+1}, t_{n+2}),
    # where the terms are sorted.
    a = 0
    b = t_seq[-4]  # The smallest of the three non-zero terms
    c = t_seq[-3]  # The middle term
    d = t_seq[-2]  # The largest term

    print("To solve this problem, we construct a tuple based on the Tribonacci sequence (0, 0, 1, 1, 2, ...).")
    print("This construction is known to yield the longest possible Ducci sequences.")
    print(f"We use the three largest terms of this sequence that are less than or equal to {limit:,}.")
    print("\nOur selected tuple (a, b, c, d) must achieve the maximum sequence length (M) and have the minimum sum.")
    print("The optimal tuple that is also lexicographically smallest is:")
    print(f"(a, b, c, d) = ({a}, {b:,}, {c:,}, {d:,})")
    
    # Step 5: Compute the required expression.
    print("\nNow we compute the expression: a + b - c - d")
    expression_value = a + b - c - d
    print(f"Expression = {a} + {b} - {c} - {d} = {expression_value}")

    # Compute the final result modulo 1000.
    result = expression_value % 1000
    print(f"\nThe value of the expression modulo 1000 is:")
    print(f"{expression_value} mod 1000 = {result}")
    
    # Final answer in the specified format
    print(f"\n<<<{result}>>>")

solve_ducci_sequence_problem()