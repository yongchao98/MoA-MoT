import math

def calculate_expected_time():
    """
    Calculates the expected number of trials to see the sequence "TENETENET"
    from a random stream of English letters.
    """
    S = "TENETENET"
    N = 26  # Size of the alphabet (English letters)
    L = len(S)

    print(f"Target sequence: S = \"{S}\" (Length L = {L})")
    print(f"Alphabet size: N = {N}")
    print("\nWe use the formula E = sum(N^k) for all k where the prefix of S of length k matches the suffix of S of length k.")
    print("\nFinding overlaps (prefix == suffix):")

    total_expectation = 0
    term_values = []
    equation_parts = []

    for k in range(1, L + 1):
        prefix = S[:k]
        suffix = S[L - k:]
        if prefix == suffix:
            print(f"  k = {k}: Prefix '{prefix}' == Suffix '{suffix}'. This is an overlap.")
            term_val = N**k
            total_expectation += term_val
            term_values.append(str(term_val))
            equation_parts.append(f"{N}^{k}")
        else:
            print(f"  k = {k}: Prefix '{prefix}' != Suffix '{suffix}'. No overlap.")
    
    # Format the final output equation as requested
    final_equation_str = " + ".join(term_values)
    final_equation_powers = " + ".join(equation_parts)

    print("\nThe expected time E is the sum of terms for each overlap found.")
    print(f"\nE = {final_equation_powers}")
    print(f"E = {final_equation_str}")
    print(f"E = {total_expectation}")

calculate_expected_time()